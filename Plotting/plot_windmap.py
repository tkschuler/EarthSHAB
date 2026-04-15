import windmap


def plot_windmap():
    wm = windmap.Windmap()
    wm.plotWind2(wm.hour_index, wm.LAT, wm.LON)