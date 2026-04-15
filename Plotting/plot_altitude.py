import matplotlib.pyplot as plt
import seaborn as sns


def plot_altitude(time_local, el, balloon_trajectory, df):
    sns.set_palette("muted")

    fig, ax = plt.subplots()
    ax.plot(time_local, el, label="reforecasted simulation")
    plt.xlabel('Datetime (MST)')
    plt.ylabel('Elevation (m)')

    if balloon_trajectory is not None:
        ax.plot(df["time"], df["altitude"], label="trajectory")

    plt.legend(loc='upper right')
    plt.title('Solar Balloon Altitude Profile')