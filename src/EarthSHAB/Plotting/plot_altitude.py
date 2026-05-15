import matplotlib as mpl
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


def plot_altitude_family(time_locals, elevations, masses, colors):
    sns.set_style("darkgrid")
    plt.rcParams.update({'font.size': 14})
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    for j, (time_local, el) in enumerate(zip(time_locals, elevations)):
        ax.plot(time_local, el, mpl.colors.rgb2hex(colors[j]), label=f"{masses[j]} kg")

    plt.xlabel('Datetime (MST)')
    plt.ylabel('Elevation (m)')
    ax.get_xaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
    ax.grid(visible=True, which='major', color='w', linewidth=1.0)
    ax.grid(visible=True, which='minor', color='w', linewidth=0.5)
    plt.legend()