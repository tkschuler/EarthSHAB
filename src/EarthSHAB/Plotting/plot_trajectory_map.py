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
    gfs,
    sim,
    t,
    start,
):
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
            gfs.LAT_LOW,
            sim.wrap_lon(gfs.lon).min(),
            gfs.LAT_HIGH,
            sim.wrap_lon(gfs.lon).max()
        )

    elif forecast_type == "ERA5":
        region = zip(*[
            (gfs.LAT_LOW, gfs.LON_LOW),
            (gfs.LAT_HIGH, gfs.LON_LOW),
            (gfs.LAT_HIGH, gfs.LON_HIGH),
            (gfs.LAT_LOW, gfs.LON_HIGH)
        ])
        gmap1.plot(lat, lon, 'red', edge_width=2.5)
        gmap1.text(coord["lat"] - .2, coord["lon"] - .2, 'Simulated Trajectory with ERA5 Reanalysis', color='red')
        gmap1.polygon(*region, color='orange', edge_width=1, alpha=.15)

    if balloon_trajectory is not None:
        if forecast_type == "GFS":
            gmap1.draw(
                "src/EarthSHAB/trajectories/" + trajectory_name + "_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html"
            )
        elif forecast_type == "ERA5":
            gmap1.draw(
                "src/EarthSHAB/trajectories/" + trajectory_name + "_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html"
            )
    else:
        if forecast_type == "GFS":
            gmap1.draw(
                "src/EarthSHAB/trajectories/PREDICTION_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html"
            )
        elif forecast_type == "ERA5":
            gmap1.draw(
                "src/EarthSHAB/trajectories/PREDICTION_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html"
            )