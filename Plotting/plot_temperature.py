import matplotlib.pyplot as plt


def plot_temperature(time_local, T_s, T_i, T_atm):
    fig, ax = plt.subplots()
    ax.plot(time_local, T_s, label="Surface Temperature")
    ax.plot(time_local, T_i, label="Internal Temperature")
    ax.plot(time_local, T_atm, label="Atmospheric Temperature")

    plt.xlabel('Datetime (MST)')
    plt.ylabel('Temperature (K)')
    plt.legend(loc='upper right')
    plt.title('Solar Balloon Temperature Profile')