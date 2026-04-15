import matplotlib.pyplot as plt
import numpy as np


def windVectorToBearing(u, v):
    bearing = np.arctan2(v, u)
    speed = np.power((np.power(u, 2) + np.power(v, 2)), 0.5)
    return [bearing, speed]


def plot_wind_bearing(time_local, x_winds_old, y_winds_old, x_winds_new, y_winds_new):
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