=====================
:mod:`Forecast`
=====================

The unified forecast reader. A single :class:`Forecast` class reads the v2
canonical netCDF schema (see the `v2 canonical schema
<https://github.com/tkschuler/EarthSHAB/blob/main/docs/forecast-schema-v2.md>`_)
regardless of whether the file originated from GFS or ERA5; the source is resolved from the
file's ``institution`` global attribute and exposed as ``Forecast.source``.

.. automodule:: Forecast
.. autoclass:: Forecast
	:members:
	:special-members: __init__
