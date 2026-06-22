from EarthSHAB.simulate import BalloonSimulation

from EarthSHAB.Plotting.plot_altitude import plot_altitude
from EarthSHAB.Plotting.plot_temperature import plot_temperature
from EarthSHAB.Plotting.plot_wind_components import plot_wind_components
from EarthSHAB.Plotting.plot_wind_bearing import plot_wind_bearing
from EarthSHAB.Plotting.plot_trajectory_map import plot_map
from EarthSHAB.Plotting.plot_windmap import plot_windmap

import matplotlib.pyplot as plt
import EarthSHAB.config as config_earth


def main():
    sim = BalloonSimulation()

    sim.run_simulation(run_reforecast=False)
    if config_earth.balloon_trajectory is not None:
        sim.run_simulation(run_reforecast=True) 
    sim_state = sim.sim_state

    # Save Variables

    scriptstartTime = sim_state["scriptstartTime"]
    coord = sim_state["coord"]
    start = sim_state["start"]
    t = sim_state["t"]
    balloon_trajectory = sim_state["balloon_trajectory"]
    forecast_type = sim_state["forecast_type"]
    trajectory_name = sim_state["trajectory_name"]
    df = sim_state["df"]

    time_local = sim_state["time_local"]
    el = sim_state["el"]
    T_s = sim_state["T_s"]
    T_i = sim_state["T_i"]
    T_atm = sim_state["T_atm"]

    x_winds_old = sim_state["x_winds_old"]
    y_winds_old = sim_state["y_winds_old"]
    x_winds_new = sim_state["x_winds_new"]
    y_winds_new = sim_state["y_winds_new"]

    lat = sim_state["lat"]
    lon = sim_state["lon"]
    lat_aprs_gps = sim_state["lat_aprs_gps"]
    lon_aprs_gps = sim_state["lon_aprs_gps"]

    gmap1 = sim_state["gmap1"]
    forecast = sim_state["forecast"]


    # Example Plots

    plot_altitude(
        time_local=time_local,
        el=el,
        balloon_trajectory=balloon_trajectory,
        df=df,
    )

    plot_temperature(
        time_local=time_local,
        T_s=T_s,
        T_i=T_i,
        T_atm=T_atm,
    )

    plot_wind_components(
        time_local=time_local,
        x_winds_new=x_winds_new,
        x_winds_old=x_winds_old,
        y_winds_new=y_winds_new,
        y_winds_old=y_winds_old,
    )

    plot_wind_bearing(
        time_local=time_local,
        x_winds_old=x_winds_old,
        y_winds_old=y_winds_old,
        x_winds_new=x_winds_new,
        y_winds_new=y_winds_new,
    )

    plot_map(
        gmap1=gmap1,
        coord=coord,
        forecast_type=forecast_type,
        balloon_trajectory=balloon_trajectory,
        trajectory_name=trajectory_name,
        lat=lat,
        lon=lon,
        lat_aprs_gps=lat_aprs_gps,
        lon_aprs_gps=lon_aprs_gps,
        df=df,
        forecast=forecast,
        sim=sim,
        t=t,
        start=start,
    )

    plot_windmap()

    plt.show()


if __name__ == "__main__":
    main()