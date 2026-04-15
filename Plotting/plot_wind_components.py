import matplotlib.pyplot as plt


def plot_wind_components(time_local, x_winds_new, x_winds_old, y_winds_new, y_winds_old):
    plt.figure()
    plt.plot(time_local, x_winds_new, label="X winds New", color="blue")
    plt.plot(time_local, x_winds_old, label="X winds Old", color="cyan")

    plt.plot(time_local, y_winds_new, label="Y winds New", color="red")
    plt.plot(time_local, y_winds_old, label="Y winds Old", color="orange")

    plt.legend(loc='upper right')
    plt.title('Wind Interpolation Comparison')