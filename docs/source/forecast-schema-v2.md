# Canonical NetCDF Forecast Schema v2

EarthSHAB v2.0 standardizes on a single netCDF forecast format: what the
**Copernicus Climate Data Store (CDS) API returns post-September 2024** for
ERA5 pressure-level reanalysis. Any new reader, downloader, or converter must
produce files that match this spec.

:::{note}
Adopting the CDS post-September-2024 layout as canonical means current ERA5
downloads load with no conversion, and GFS is the only source that needs a
converter (`saveNETCDF.py` for live forecasts, `saveNETCDF_archive.py` for
historical cycles from the AWS GFS archive — both emit identical v2 output).

For converting between older/alternate forecast layouts (e.g. pre-September-2024
ERA5 or other variable/dimension naming), see the
[ERA5-Utils toolkit](https://github.com/tkschuler/ERA5-Utils).
:::

Reference files (the bundled SHAB14-V flight):
`src/EarthSHAB/forecasts/SHAB14V_ERA5_20220822_20220823.nc` (ERA5) and
`src/EarthSHAB/forecasts/gfs_0p25_20220822_12.nc` (GFS, migrated from v1 to v2 canonical format).

---

## 1. Global attributes

| Attribute       | Required | Value                                                                                  |
|-----------------|----------|----------------------------------------------------------------------------------------|
| `Conventions`   | yes      | `"CF-1.7"` (string)                                                                    |
| `institution`   | yes      | `"NOAA/NCEP (GFS)"` for GFS files; `"ECMWF (ERA5)"` for ERA5 files                     |
| `history`       | yes      | One line of provenance (timestamp + tool that wrote the file)                          |
| `GRIB_centre`   | optional | Present on direct Copernicus downloads (`"ecmf"`). Tolerated; not required             |

:::{warning}
The reader's format-detection uses `Conventions` + variable names. Any file
lacking `Conventions == "CF-1.7"` AND containing the old GFS variables
(`ugrdprs`, `vgrdprs`, `hgtprs`) is treated as a v1 file and refused with a
migration message (see [migration-v2.md](migration-v2.md)).
:::

### Automatic Netcdf Forecast Sourcing 

`Forecast.source` (the field used for plot labels and evaluation grouping) is
resolved from `institution`:

* contains `"NOAA"`, `"GFS"`, or `"NCEP"` → `"GFS"`
* contains `"ECMWF"` or `"ERA5"`         → `"ERA5"`
* missing/empty                          → fall back to filename pattern
  (basename contains `"gfs"` → `"GFS"`; contains `"era5"` → `"ERA5"`;
  otherwise `"unknown"`)

`saveNETCDF.py`, `saveNETCDF_archive.py`, and `migrate_v1.py` all write the appropriate `institution`
value. Files migrated before this convention was introduced carry an empty
`institution` and rely on the filename fallback.

### Storage convention

**Every v2 file is a bounded subset.** No full-world arrays with
mask-based subsetting; the shape of `u/v/z/t` IS the data extent. Missing
samples inside the subset (e.g., at high pressure levels above the model top)
are represented as `_FillValue` / NaN and resolved by the reader's
`fill_missing_data` 1-D interpolation. The reader does not perform any
mask-based outer-bounding-box detection.

This is enforced upstream: `saveNETCDF.py` downloads via NOAA's GRIB-filter
`subregion`, `saveNETCDF_archive.py` crops each archive step to the configured
bounding box, `migrate_v1.py` strips masked padding from v1 GFS archives, and
raw Copernicus CDS ERA5 downloads are already subsets.

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

:::{note}
GFS downloads natively use 0-360; `saveNETCDF.py` must convert to -180-180 before writing.
:::

---

## 4. Scalar / auxiliary coordinates

Present on real CDS files; the canonical schema TOLERATES them but the
reader treats them as informational only.

| Variable | Shape           | Purpose                                                  |
|----------|-----------------|----------------------------------------------------------|
| `number` | scalar (int64)  | Ensemble member ID; always `0` for our single-realization use case |
| `expver` | `(valid_time,)` string | ECMWF experiment version (`"0001"` for operational ERA5) |

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

:::{important}
`z` is **geopotential** (m²/s²), not geopotential height (m). The reader does
the `/g` conversion when interpolating against an altitude query.
:::

---

## 6. Index axis summary

For a query at `(t_query, alt_query_m, lat_query, lon_query)` the reader:

1. finds nearest indices: `valid_time` (ascending), `latitude` (descending),
   `longitude` (ascending; query in -180-180);
2. builds the altitude column `z[t_idx, :, lat_idx, lon_idx] / g`;
3. interpolates `u`, `v` over altitude (per `wind_interpolation`), then over
   `valid_time` between the two enclosing time steps.

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

`saveNETCDF.py` (current) and `saveNETCDF_archive.py` (historical, AWS archive) both do all the conversion at download time; pre-v2 archived files are migrated by the `migrate_v1` CLI.

:::{seealso}
[migration-v2.md](migration-v2.md) — detecting v1 files, the converter, and rollback.
:::

---

## 9. Reference

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