import EarthSHAB.windmap as windmap


def plot_windmap(interpolation='linear'):
    """Plot a 3D windrose.

    :param interpolation: ``'linear'`` interpolates speed/direction between adjacent pressure
        levels with angle-wrap correction (default); ``'spline'`` fits a CubicSpline across all
        pressure levels for a smoother continuous curve
    :type interpolation: str
    """
    wm = windmap.Windmap()
    wm.plotWind2(wm.hour_index, wm.LAT, wm.LON, interpolation=interpolation)