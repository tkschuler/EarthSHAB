# Migrating archived forecasts to v2

EarthSHAB v2.0 unifies the GFS and ERA5 readers around a single canonical
netCDF layout (see [forecast-schema-v2.md](forecast-schema-v2.md)). Files in
the older v1 formats are **rejected on load** with `ForecastFormatError`.
Convert them once with the migrate_v1 CLI; the original is preserved as a
`.v1.nc` sibling for rollback.

## TL;DR

```bash
python -m EarthSHAB.forecast_processing.migrate_v1 src/EarthSHAB/forecasts/
```

Per file, the script:

1. Renames `<name>.nc` → `<name>.v1.nc` (backup).
2. Writes a canonical-format file in its place at `<name>.nc`.

The operation is idempotent. Re-running on the same directory only emits
`SKIP` lines. A `--dry-run` flag is available.

## What the reader does when it sees a v1 file

```
ForecastFormatError: forecasts/foo.nc: detected v1 GFS (ugrdprs/vgrdprs/hgtprs
variables); EarthSHAB v2.0+ requires the canonical schema documented in
docs/forecast-schema-v2.md. Convert this file in place with:
    python -m EarthSHAB.forecast_processing.migrate_v1 forecasts/foo.nc
The original will be backed up alongside it as <name>.v1.nc. See
docs/migration-v2.md for details.
```

The reader prints the exact command to run, with the offending path baked in.

## Which legacy formats are recognised

| Tag       | Signature                                                        | Source                                                              |
|-----------|------------------------------------------------------------------|---------------------------------------------------------------------|
| `v1_gfs`  | Variables `ugrdprs`, `vgrdprs`, `hgtprs`, `tmpprs`               | Output of pre-v2 `saveNETCDF.py`                                    |
| `v1_era5` | Vars `u`,`v`,`z`,`t` with dims `time`,`level`,`latitude`,`longitude` | Output of the removed `convert_era52gfs.py`; pre-Sep-2024 raw CDS downloads |

Anything else exits with `ERROR`. v1 ERA5 covers both the CF-1.7 "processed"
intermediate and CF-1.6 pre-Sep-2024 CDS downloads — xarray decodes either
time encoding transparently, and both layouts need the same axis fixes.

## What the converter changes

**v1 GFS → v2 canonical**

- Renames dims `time`→`valid_time`, `lev`→`pressure_level`, `lat`→`latitude`, `lon`→`longitude`.
- Renames variables `ugrdprs`→`u`, `vgrdprs`→`v`, `tmpprs`→`t`, `hgtprs`→`z`.
- Converts geopotential height `hgtprs` (gpm) to geopotential `z = h * g` (m²/s²).
- Converts longitude from `[0, 360)` to `[-180, 180)` and re-sorts.
- Flips latitude from ascending to descending.
- Re-encodes time as `seconds since 1970-01-01` / `proleptic_gregorian`. The
  reader applies the same `has_year_zero=True` decode the legacy `GFS.py` used,
  so the resulting wall-clock timestamps match what the v1 reader surfaced.

**v1 ERA5 → v2 canonical**

- Renames dims `time`→`valid_time`, `level`→`pressure_level`.
- Flips `level` from ascending hPa to descending hPa (data reverses with it).
- Re-encodes time as `seconds since 1970-01-01` / `proleptic_gregorian` (a
  no-op on already-canonical encodings; a real conversion on the CF-1.6
  "hours since 1900-01-01" pre-Sep-2024 downloads).
- Latitude order and longitude convention are already correct; `z` is already
  geopotential. Both are left as-is.

All canonical files are written with `Conventions = "CF-1.7"`, float32 data,
and the dimension order `(valid_time, pressure_level, latitude, longitude)`.

## Why no auto-conversion in the reader

Silent on-load conversion would mutate user data without an audit trail, hide
the cost of the change, and make it impossible to reproduce a prior run. The
reader fails loudly with the exact CLI to run; the converter writes a `.v1.nc`
backup so the original is always recoverable.

## CLI reference

```
python -m EarthSHAB.forecast_processing.migrate_v1 [--dry-run] [--verbose] <target> [<target> ...]
```

- `<target>` is a `.nc` file or a directory. Directories are scanned
  non-recursively for `*.nc` files; existing `*.v1.nc` backups are skipped.
- Default output is a `tqdm` progress bar plus a final summary line; only
  `ERROR` lines stream to stdout. Pass `--verbose` (`-v`) to print one
  status line per file (`OK` / `SKIP` / `ERROR` / `DRY-RUN`).
- Exit status is `0` if every target was OK or SKIP, `1` if any ERROR.
- Status prefixes are stable for tooling: `OK`, `SKIP`, `ERROR`, `DRY-RUN`.

## Verifying the migration

After running the converter, confirm no v1-format files remain:

```bash
python -m EarthSHAB.forecast_processing.verify_migration src/EarthSHAB/forecasts/
```

The verifier walks the same targets, inspects each `.nc` with the same
`detect_format` logic the migrator uses, and reports any files that are still
in `v1_gfs` or `v1_era5` layout or are structurally `unknown`. `.v1.nc`
backups are skipped (counted separately). Exit status is `0` if every
inspected file is `v2`, `1` if any file is v1 / unknown / can't be opened.

```
python -m EarthSHAB.forecast_processing.verify_migration [--verbose] <target> [<target> ...]
```

- Default output is a `tqdm` progress bar plus a final summary line; only
  non-`v2` findings stream to stdout. Pass `--verbose` (`-v`) to print one
  line per file.
- A clean run ends with `OK: all inspected forecast files are v2 canonical.`
- A failed run ends with `FAIL: N file(s) need attention. Re-run migrate_v1
  on the listed paths.`

## Rolling back

Delete the converted file and rename the backup:

```bash
rm forecasts/foo.nc
mv forecasts/foo.v1.nc forecasts/foo.nc
```

Subsequent loads will fail with `ForecastFormatError` again, which is the
correct signal that the file is back in the unsupported v1 state.
