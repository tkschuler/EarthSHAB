# Canonical Forecast Schema (v2)

EarthSHAB v2.0 standardizes on a single netCDF forecast format, matching
what the **Copernicus Climate Data Store (CDS) API returns post-September
2024** when downloading ERA5 pressure-level reanalysis. The reference file
reference files are the bundled SHAB14-V v2 forecasts
`src/EarthSHAB/forecasts/SHAB14V_ERA5_20220822_20220823.nc` (ERA5) and
`src/EarthSHAB/forecasts/gfs_0p25_20220822_12.nc` (GFS, migrated to canonical) —
any new reader, downloader, or converter must produce files that match this
spec.

> **Why this format?** It is what the upstream Copernicus API now returns by
> default. Adopting it as the canonical schema means archived ERA5 downloads
> work without conversion, future ERA5 downloads work as-is, and GFS becomes
> the only data source that requires a converter (saveNETCDF.py).

---

## 1. Global attributes

| Attribute       | Required | Value                                                                                  |
|-----------------|----------|----------------------------------------------------------------------------------------|
| `Conventions`   | yes      | `"CF-1.7"` (string)                                                                    |
| `institution`   | yes      | `"NOAA/NCEP (GFS)"` for GFS files; `"ECMWF (ERA5)"` for ERA5 files                     |
| `history`       | yes      | One line of provenance (timestamp + tool that wrote the file)                          |
| `GRIB_centre`   | optional | Present on direct Copernicus downloads (`"ecmf"`). Tolerated; not required             |

The reader's format-detection logic uses `Conventions` + variable names. Any
file lacking `Conventions == "CF-1.7"` AND containing the old GFS variables
(`ugrdprs`, `vgrdprs`, `hgtprs`) is treated as a v1 file and refused with a
migration message.

### Source provenance

`Forecast.source` (the field used for plot labels and evaluation grouping) is
resolved from `institution`:

* contains `"NOAA"`, `"GFS"`, or `"NCEP"` → `"GFS"`
* contains `"ECMWF"` or `"ERA5"`         → `"ERA5"`
* missing/empty                          → fall back to filename pattern
  (basename contains `"gfs"` → `"GFS"`; contains `"era5"` → `"ERA5"`;
  otherwise `"unknown"`)

`saveNETCDF.py` and `migrate_v1.py` both write the appropriate `institution`
value. Files already migrated prior to this convention (Phase 5) carry an
empty `institution` and rely on the filename fallback.

### Storage convention

**Every v2 file is a tight bounding-box subset.** No full-world arrays with
mask-based subsetting; the shape of `u/v/z/t` IS the data extent. Missing
samples inside the subset (e.g., at high pressure levels above the model top)
are represented as `_FillValue` / NaN and resolved by the reader's
`fill_missing_data` 1-D interpolation. The reader does not perform any
mask-based outer-bounding-box detection.

This convention is enforced upstream:

* `saveNETCDF.py` downloads via NOAA's GRIB filter `subregion` parameter, so
  freshly downloaded GFS files are already subsets.
* `migrate_v1.py`'s `_valid_subset_indices()` strips the masked padding from
  v1 GFS global archives during conversion.
* Raw Copernicus CDS ERA5 downloads have always been subsets.

---

## 2. Dimensions

Four dimensions, in this exact order on every data variable:

| Dimension        | Type      | Notes                                              |
|------------------|-----------|----------------------------------------------------|
| `valid_time`     | int64     | UTC instants. **Renamed from v1 `time`**.          |
| `pressure_level` | float64   | Isobaric levels in hPa. **Renamed from v1 `level`**. |
| `latitude`       | float64   | Degrees north.                                     |
| `longitude`      | float64   | Degrees east.                                      |

Dimension ORDER on data vars: `(valid_time, pressure_level, latitude, longitude)`.

---

## 3. Coordinate variables

### 3.1 `valid_time` (1D, dim `valid_time`)

| Attribute       | Value                                  |
|-----------------|----------------------------------------|
| `units`         | `"seconds since 1970-01-01"`           |
| `calendar`      | `"proleptic_gregorian"`                |
| `standard_name` | `"time"`                               |
| `long_name`     | `"time"`                               |
| dtype           | `int64`                                |
| stored order    | strictly ascending                     |

> **No `has_year_zero` quirk.** Using `proleptic_gregorian` with epoch
> seconds round-trips cleanly to stdlib `datetime`; no 2-day calendar offset.

### 3.2 `pressure_level` (1D, dim `pressure_level`)

| Attribute          | Value                |
|--------------------|----------------------|
| `units`            | `"hPa"`              |
| `positive`         | `"down"`             |
| `stored_direction` | `"decreasing"`       |
| `standard_name`    | `"air_pressure"`     |
| `long_name`        | `"pressure"`         |
| `_FillValue`       | `NaN`                |
| dtype              | `float64`            |
| stored order       | **DESCENDING** in hPa (e.g., `[1000, 975, …, 2, 1]`) → ASCENDING in altitude |

Real CDS files publish 37 standard pressure levels: `[1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100, 70, 50, 30, 20, 10, 7, 5, 3, 2, 1]`. Subsets are allowed; ordering must still be descending hPa.

### 3.3 `latitude` (1D, dim `latitude`)

| Attribute          | Value                |
|--------------------|----------------------|
| `units`            | `"degrees_north"`    |
| `standard_name`    | `"latitude"`         |
| `long_name`        | `"latitude"`         |
| `stored_direction` | `"decreasing"`       |
| `_FillValue`       | `NaN`                |
| dtype              | `float64`            |
| stored order       | **DESCENDING** (e.g., `[51.0, 50.75, …, 21.25]`) — north to south |
| grid spacing       | uniform; 0.25° for ERA5; 0.25° for GFS at 0p25 resolution            |

### 3.4 `longitude` (1D, dim `longitude`)

| Attribute       | Value                |
|-----------------|----------------------|
| `units`         | `"degrees_east"`     |
| `standard_name` | `"longitude"`        |
| `long_name`     | `"longitude"`        |
| `_FillValue`    | `NaN`                |
| dtype           | `float64`            |
| stored order    | **ASCENDING** (e.g., `[-121.75, -121.5, …, -72.0]`) |
| convention      | **`-180` to `180`** (NOT `0` to `360`)              |

> GFS downloads natively use 0-360; saveNETCDF.py must convert to -180 to 180 before writing.

---

## 4. Scalar / auxiliary coordinates

Present on real CDS files; the canonical schema TOLERATES them but the
reader treats them as informational only.

| Variable | Shape           | Purpose                                                  |
|----------|-----------------|----------------------------------------------------------|
| `number` | scalar (int64)  | Ensemble member ID; always `0` for our single-realization use case |
| `expver` | `(valid_time,)` string | ECMWF experiment version (`"0001"` for operational ERA5) |

The reader MUST NOT depend on these being present. The GFS converter is NOT required to emit them.

---

## 5. Data variables

Four required data variables. All have:
- **dtype**: `float32`
- **dims**: `(valid_time, pressure_level, latitude, longitude)`
- **`_FillValue`**: `NaN` (float32)
- **`coordinates`** attribute listing applicable scalars: `"number valid_time isobaricInhPa latitude longitude expver"` is the CDS form. The reader accepts any value (or missing).

| Var | `standard_name`     | `long_name`            | `units`     | Notes                                              |
|-----|---------------------|------------------------|-------------|----------------------------------------------------|
| `u` | `eastward_wind`     | `U component of wind`  | `m s**-1`   | Eastward (positive = wind blowing TO the east)     |
| `v` | `northward_wind`    | `V component of wind`  | `m s**-1`   | Northward (positive = wind blowing TO the north)   |
| `z` | `geopotential`      | `Geopotential`         | `m**2 s**-2`| **NOT altitude in meters.** Divide by `g = 9.80665` to obtain geopotential height in meters |
| `t` | `air_temperature`   | `Temperature`          | `K`         | Optional in v1 ERA5, required in v2                |

`z` is **geopotential** (m²/s²), not geopotential height (m). The reader is responsible for the `/g` conversion when interpolating against an altitude query.

GRIB attribute pass-throughs (`GRIB_paramId`, `GRIB_centre`, …) are TOLERATED but not required. The GFS converter need not emit them.

---

## 6. Index axis summary

For a query at `(t_query, alt_query_m, lat_query, lon_query)`, the canonical reader:

1. `valid_time_idx = nearest_idx(valid_time, t_query)` (ascending search)
2. `latitude_idx  = nearest_idx(latitude, lat_query)` (descending search)
3. `longitude_idx = nearest_idx(longitude, lon_query_in_minus180_to_180)` (ascending)
4. Build column `alt_col_m = z[t_idx, :, lat_idx, lon_idx] / g`
5. Reverse the column to ascending altitude (since `pressure_level` is descending hPa → altitude ascends with index? **No** — descending hPa means highest pressure first = LOWEST altitude first, so the column IS ASCENDING in altitude as level index grows). **No reversal needed.**
6. Interpolate `u`, `v` over altitude (per the selected `wind_interpolation` method)
7. Interpolate over `valid_time` between the two enclosing time indices

> Sanity check: real ERA5 `z[0, :, 0, 0]` decreases from ~475809 m²/s² (level=1 hPa, ~48km) to ~surface as level index increases from 0. Because `pressure_level` is DESCENDING hPa with index, the FIRST level index (= 1000 hPa) is the LOWEST altitude. So as level index grows, pressure DECREASES and altitude INCREASES. The altitude column IS ascending with index, so the `Forecast` reader does no reversal. (The legacy `ERA5.py` applied an unnecessary `[::-1]` reversal — a v1-format quirk that was removed when the readers collapsed into `Forecast`.)

---

## 7. Differences from v1 ("processed" ERA5 format the legacy ERA5.py expected)

| Item                | v1 (legacy ERA5.py)                   | v2 (canonical)                            |
|---------------------|---------------------------------------|-------------------------------------------|
| Time dim/var name   | `time`                                | `valid_time`                              |
| Level dim/var name  | `level`                               | `pressure_level`                          |
| Time units          | varies                                | `seconds since 1970-01-01`                |
| Time calendar       | varies                                | `proleptic_gregorian`                     |
| Pressure ordering   | ascending hPa                         | **descending hPa** (level idx 0 = 1000 hPa) |
| `z` semantics       | already-divided-by-g (meters)          | **geopotential (m²/s²)**; reader divides  |
| `t` (temperature)   | optional                              | required                                  |
| CF Conventions attr | absent or older                       | `CF-1.7`                                  |

The Forecast class auto-detects v1 vs v2 and refuses v1 with a clear migration message (see [migration-v2.md](migration-v2.md)).

---

## 8. Differences from v1 GFS format (saveNETCDF.py output)

| Item                       | v1 GFS (`gfs_0p25_*.nc`)                                  | v2 (canonical)                         |
|----------------------------|------------------------------------------------------------|----------------------------------------|
| Variable names             | `ugrdprs`, `vgrdprs`, `hgtprs`, `tmpprs`                  | `u`, `v`, `z`, `t`                      |
| Dimension names            | `time`, `lev`, `lat`, `lon`                                | `valid_time`, `pressure_level`, `latitude`, `longitude` |
| Time encoding              | `days since 0001-01-01`, no units attr (hardcoded in reader) | `seconds since 1970-01-01`, calendar attr |
| Altitude variable          | `hgtprs` (meters)                                          | `z` (m²/s², geopotential)              |
| Longitude convention       | `[0, 360)`                                                 | `[-180, 180)`                          |
| Latitude order             | ascending                                                  | descending                             |
| Level (`lev`) order        | descending hPa                                             | same (`pressure_level` descending hPa) |
| `Conventions` attr         | absent                                                     | `CF-1.7`                               |

`saveNETCDF.py` is the place to do all the conversion at download time; archived files are migrated by the `migrate_v1` CLI (see Changelog v2.0.0 Phase 3).

---

## 9. Validation

The canonical Forecast class validates on `__init__` and raises `ForecastFormatError` with a specific message if any of these fail:

1. Required dimensions (`valid_time`, `pressure_level`, `latitude`, `longitude`) all present with `len >= 1`.
2. Required data variables (`u`, `v`, `z`, `t`) all present with `dtype == float32` and correct dim ordering.
3. `pressure_level` strictly descending; `valid_time` strictly ascending; `latitude` strictly descending; `longitude` strictly ascending.
4. `valid_time.units` parseable as a CF time unit; `valid_time.calendar` set.
5. `Conventions == "CF-1.7"` (warn if missing or different value but still attempt to load — interop with future CF revisions).
6. `longitude` values all within `[-180, 180]`.

A separate `migrate_v1` script handles detection of the two v1 formats (by variable-name signature) and produces a canonical file in-place with the original backed up to `<name>.v1.nc`.

---

## 10. Reference

The canonical example files are the bundled SHAB14-V v2 forecasts:

```
src/EarthSHAB/forecasts/SHAB14V_ERA5_20220822_20220823.nc   # ERA5
src/EarthSHAB/forecasts/gfs_0p25_20220822_12.nc             # GFS (migrated)
```

Inspect either directly with:

```python
import netCDF4
ds = netCDF4.Dataset('src/EarthSHAB/forecasts/SHAB14V_ERA5_20220822_20220823.nc')
print(ds)
```

When in doubt about an edge case not covered above, the reference files'
behavior is authoritative.