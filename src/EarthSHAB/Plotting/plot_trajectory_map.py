import os

if not os.path.exists('src/EarthSHAB/trajectories'):
    os.makedirs('src/EarthSHAB/trajectories')


def draw_bounding_box(gmap1, min_lat, min_lon, max_lat, max_lon):
    lats = []
    lons = []

    # Bottom edge
    for lon in range(int(min_lon), int(max_lon) + 1, 10):
        lats.append(min_lat)
        lons.append(lon)

    # Right edge
    for lat in range(int(min_lat), int(max_lat) + 1, 10):
        lats.append(lat)
        lons.append(max_lon)

    # Top edge
    for lon in range(int(max_lon), int(min_lon) - 1, -10):
        lats.append(max_lat)
        lons.append(lon)

    # Left edge
    for lat in range(int(max_lat), int(min_lat) - 1, -10):
        lats.append(lat)
        lons.append(min_lon)

    gmap1.polygon(lats, lons, 'cornflowerblue', edge_width=5, alpha=.2)


def _drop_nan_coords(lats, lons):
    """Return (lats, lons) with any pair containing NaN removed."""
    import math
    def _isnan(v):
        try:
            return math.isnan(v)
        except (TypeError, ValueError):
            return False
    pairs = [(la, lo) for la, lo in zip(lats, lons) if not _isnan(la) and not _isnan(lo)]
    if not pairs:
        return [], []
    return zip(*pairs)


def plot_map(
    gmap1,
    coord,
    forecast_type,
    balloon_trajectory,
    trajectory_name,
    lat,
    lon,
    lat_aprs_gps,
    lon_aprs_gps,
    df,
    forecast,
    sim,
    t,
    start,
    html_prefix=None,
    output_dir="src/EarthSHAB/trajectories/",
):
    lat_aprs_gps, lon_aprs_gps = _drop_nan_coords(lat_aprs_gps, lon_aprs_gps)
    lat, lon = _drop_nan_coords(lat, lon)

    if balloon_trajectory is not None:
        if forecast_type == "GFS":
            gmap1.plot(lat_aprs_gps, lon_aprs_gps, 'cyan', edge_width=2.5)
            gmap1.text(
                coord["lat"] - .3,
                coord["lon"] - .2,
                trajectory_name + " Alt + " + forecast_type + " Wind Data",
                color='cyan'
            )
        elif forecast_type == "ERA5":
            gmap1.plot(lat_aprs_gps, lon_aprs_gps, 'orange', edge_width=2.5)
            gmap1.text(
                coord["lat"] - .3,
                coord["lon"] - .2,
                trajectory_name + " Alt + " + forecast_type + " Wind Data",
                color='orange'
            )

    if forecast_type == "GFS":
        gmap1.plot(lat, lon, 'blue', edge_width=2.5)
        gmap1.text(coord["lat"] - .2, coord["lon"] - .2, 'Simulated Trajectory with GFS Forecast', color='blue')
        draw_bounding_box(
            gmap1,
            forecast.LAT_LOW,
            sim.wrap_lon(forecast.lon).min(),
            forecast.LAT_HIGH,
            sim.wrap_lon(forecast.lon).max()
        )

    elif forecast_type == "ERA5":
        region = zip(*[
            (forecast.LAT_LOW, forecast.LON_LOW),
            (forecast.LAT_HIGH, forecast.LON_LOW),
            (forecast.LAT_HIGH, forecast.LON_HIGH),
            (forecast.LAT_LOW, forecast.LON_HIGH)
        ])
        gmap1.plot(lat, lon, 'red', edge_width=2.5)
        gmap1.text(coord["lat"] - .2, coord["lon"] - .2, 'Simulated Trajectory with ERA5 Reanalysis', color='red')
        gmap1.polygon(*region, color='orange', edge_width=1, alpha=.15)

    prefix = (html_prefix or "") + (trajectory_name if balloon_trajectory is not None else "PREDICTION")

    # Encode the forecast's spatial grid and temporal cadence in the filename,
    # e.g. "GFS_0p25deg_3hr". Spatial resolution comes from the grid spacing,
    # temporal step from the file cadence (Forecast.resolution_hr).
    res_token = ""
    if len(forecast.lat) >= 2:
        res_deg = abs(float(forecast.lat[1]) - float(forecast.lat[0]))
        res_token = f"_{res_deg:g}deg".replace(".", "p")
    step_hr = getattr(forecast, "resolution_hr", None)
    step_token = f"_{step_hr:g}hr".replace(".", "p") if step_hr is not None else ""

    os.makedirs(output_dir, exist_ok=True)
    gmap1.draw(
        os.path.join(
            output_dir,
            f"{prefix}_{forecast_type}{res_token}{step_token}_{t.year}_{t.month}_{start.day}.html"
        )
    )