from termcolor import colored
import matplotlib.pyplot as plt
import time as tm
import os
import numpy as np
import seaborn as sns
import windmap


if not os.path.exists('trajectories'):
    os.makedirs('trajectories')


def windVectorToBearing(u, v):
    bearing = np.arctan2(v, u)
    speed = np.power((np.power(u, 2) + np.power(v, 2)), 0.5)
    return [bearing, speed]


def draw_bounding_box(gmap1, min_lat, min_lon, max_lat, max_lon):
    """
    Draws a rectangular bounding box on a Google Map and saves it as an HTML file.

    Args:
        min_lat (float): Minimum latitude (southern boundary).
        min_lon (float): Minimum longitude (western boundary).
        max_lat (float): Maximum latitude (northern boundary).
        max_lon (float): Maximum longitude (eastern boundary).
    """

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


def run_all_plots(sim_state, sim):
    scriptstartTime = sim_state["scriptstartTime"]

    coord = sim_state["coord"]
    start = sim_state["start"]
    t = sim_state["t"]
    balloon_trajectory = sim_state["balloon_trajectory"]
    forecast_type = sim_state["forecast_type"]
    trajectory_name = sim_state["trajectory_name"]
    df = sim_state["df"]

    T_s = sim_state["T_s"]
    T_i = sim_state["T_i"]
    T_atm = sim_state["T_atm"]
    el = sim_state["el"]

    x_winds_old = sim_state["x_winds_old"]
    y_winds_old = sim_state["y_winds_old"]
    x_winds_new = sim_state["x_winds_new"]
    y_winds_new = sim_state["y_winds_new"]

    lat = sim_state["lat"]
    lon = sim_state["lon"]
    lat_aprs_gps = sim_state["lat_aprs_gps"]
    lon_aprs_gps = sim_state["lon_aprs_gps"]

    gmap1 = sim_state["gmap1"]
    gfs = sim_state["gfs"]

    # support either new or old naming
    time_local = sim_state["time_local"]

    sns.set_palette("muted")

    fig, ax = plt.subplots()
    ax.plot(time_local, el, label="Reforecasted simulation")
    plt.title('Balloon Altitude Profile')
    plt.xlabel('Datetime (MST)')
    plt.ylabel('Elevation (m)')

    if balloon_trajectory is not None:
        ax.plot(df["time"], df["altitude"], label="trajectory")

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

    plt.legend(loc='upper right')

    fig2, ax2 = plt.subplots()
    ax2.plot(time_local, T_s, label="Surface Temperature")
    ax2.plot(time_local, T_i, label="Internal Temperature")
    ax2.plot(time_local, T_atm, label="Atmospheric Temperature")
    plt.xlabel('Datetime (MST)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc='upper right')
    plt.title('Solar Balloon Temperature Profile')

    plt.figure()
    plt.plot(time_local, x_winds_new, label="X winds New", color="blue")
    plt.plot(time_local, x_winds_old, label="X winds Old", color="cyan")
    plt.plot(time_local, y_winds_new, label="Y winds New", color="red")
    plt.plot(time_local, y_winds_old, label="Y winds Old", color="orange")
    plt.legend(loc='upper right')
    plt.title('Wind Interpolation Comparison')

    plt.figure()
    if any(x_winds_old):
        plt.plot(
            time_local,
            np.degrees(windVectorToBearing(x_winds_old, y_winds_old)[0]),
            label="Bearing old",
            color="blue"
        )
    plt.plot(
        time_local,
        np.degrees(windVectorToBearing(x_winds_new, y_winds_new)[0]),
        label="Bearing New",
        color="red"
    )
    plt.legend(loc='upper right')
    plt.title('Wind Interpolation Comparison')

    # Outline Downloaded NOAA forecast subset:
    if forecast_type == "GFS":
        gmap1.plot(lat, lon, 'blue', edge_width=2.5)
        gmap1.text(coord["lat"] - .2, coord["lon"] - .2, 'Simulated Trajectory with GFS Forecast', color='blue')
        draw_bounding_box(gmap1, gfs.LAT_LOW, sim.wrap_lon(gfs.lon).min(), gfs.LAT_HIGH, sim.wrap_lon(gfs.lon).max())

    elif forecast_type == "ERA5":
        region = zip(*[
            (gfs.LAT_LOW, gfs.LON_LOW),
            (80, gfs.LON_LOW),
            (80, gfs.LON_HIGH),
            (gfs.LAT_LOW, gfs.LON_HIGH)
        ])
        gmap1.plot(lat, lon, 'red', edge_width=2.5)
        gmap1.text(coord["lat"] - .2, coord["lon"] - .2, 'Simulated Trajectory with ERA5 Reanalysis', color='red')
        gmap1.polygon(*region, color='orange', edge_width=1, alpha=.15)

    year = str(tm.localtime()[0])
    month = str(tm.localtime()[1]).zfill(2)
    day = str(tm.localtime()[2]).zfill(2)

    if balloon_trajectory is not None:
        if forecast_type == "GFS":
            gmap1.draw("trajectories/" + trajectory_name + "_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html")
        elif forecast_type == "ERA5":
            gmap1.draw("trajectories/" + trajectory_name + "_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html")
    else:
        if forecast_type == "GFS":
            gmap1.draw("trajectories/PREDICTION_GFS_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html")
        elif forecast_type == "ERA5":
            gmap1.draw("trajectories/PREDICTION_ERA5_" + str(t.year) + "_" + str(t.month) + "_" + str(start.day) + ".html")

    executionTime = (tm.time() - scriptstartTime)
    print('\nSimulation executed in ' + str(executionTime) + ' seconds.')

    wm = windmap.Windmap()
    wm.plotWind2(wm.hour_index, wm.LAT, wm.LON)

    plt.show()