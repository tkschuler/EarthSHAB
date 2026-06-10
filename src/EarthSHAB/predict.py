"""predict.py creates a family of predictions at different float altitudes.

Float altitudes are adjusted by changing the payload mass in .25 kg increments.
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import gmplot
from termcolor import colored

import EarthSHAB.config_earth as config_earth
from EarthSHAB.simulate import BalloonSimulation
from EarthSHAB.Plotting.plot_altitude import plot_altitude_family
from EarthSHAB.Plotting.plot_trajectory_map import build_trajectory_filename
from EarthSHAB.Plotting.plot_windmap import plot_windmap


def predict():
    coord = config_earth.simulation['start_coord']
    masses = [0, .25, .5, .75, 1, 1.25, 1.5, 1.75, 2]

    cmap = mpl.colormaps['rainbow_r']
    colors = cmap(np.linspace(0, 1, len(masses)))

    gmap1 = gmplot.GoogleMapPlotter(coord["lat"], coord["lon"], 10)

    time_locals = []
    elevations = []

    for j, mass in enumerate(masses):
        print(colored("---------------------------" + str(mass) + "kg-------------------------------", "magenta"))

        config_earth.balloon_properties['mp'] = mass

        sim = BalloonSimulation()
        sim.run_simulation(run_reforecast=False, descent_correction=True)
        sim_state = sim.sim_state

        time_locals.append(sim_state["time_local"])
        elevations.append(sim_state["el"])

        gmap1.plot(sim_state["lat"], sim_state["lon"], mpl.colors.rgb2hex(colors[j]), edge_width=2.5)

    t = sim_state["t"]
    start = sim_state["start"]
    forecast = sim_state["forecast"]
    forecast_type = sim_state["forecast_type"]

    plot_altitude_family(time_locals, elevations, masses, colors)

    region = zip(*[
        (forecast.LAT_LOW, forecast.LON_LOW),
        (forecast.LAT_HIGH, forecast.LON_LOW),
        (forecast.LAT_HIGH, forecast.LON_HIGH),
        (forecast.LAT_LOW, forecast.LON_HIGH)
    ])
    gmap1.polygon(*region, color='cornflowerblue', edge_width=1, alpha=.2)

    # Match the trajectory HTML naming used by main.py / plot_map.
    output_dir = "src/EarthSHAB/trajectories/"
    os.makedirs(output_dir, exist_ok=True)
    gmap1.draw(
        os.path.join(
            output_dir,
            build_trajectory_filename(forecast, forecast_type, t, start, prefix="MULTI_PREDICTION")
        )
    )

    plot_windmap()
    plt.show()


if __name__ == "__main__":
    predict()
