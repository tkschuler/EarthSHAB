=====================
:mod:`Forecast`
=====================

The unified forecast reader. A single :class:`Forecast` class reads the v2
canonical netCDF schema (see :doc:`../forecast-schema-v2`)
regardless of whether the file originated from GFS or ERA5; the source is resolved from the
file's ``institution`` global attribute and exposed as ``Forecast.source``.

:class:`Forecast` is a thin **orchestrator**: it owns the public API —
:func:`~Forecast.Forecast.getNewCoord` (scalar) and
:func:`~Forecast.Forecast.getNewCoords` (batch) — and the trajectory advection,
while delegating the wind lookup to pluggable backends
(:mod:`EarthSHAB.utils.forecast_backends`), the altitude/time interpolation to
:mod:`EarthSHAB.utils.wind_interp`, and the advection step to
:mod:`EarthSHAB.utils.advection`. The wind field is selected by three orthogonal
config fields — ``wind_interpolation``, ``advection``, and ``backend`` — all
documented on :doc:`wind_interpolation`. These ``utils`` modules are shared
byte-for-byte with the HAB-COM project so interpolation/physics edits stay
portable between the two codebases.

.. automodule:: Forecast
.. autoclass:: Forecast
	:members:
	:special-members: __init__
