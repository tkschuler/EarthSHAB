=====================
:mod:`Forecast`
=====================

The unified forecast reader. A single :class:`Forecast` class reads the v2
canonical netCDF schema (see :doc:`../forecast-schema-v2`)
regardless of whether the file originated from GFS or ERA5; the source is resolved from the
file's ``institution`` global attribute and exposed as ``Forecast.source``.

.. automodule:: Forecast
.. autoclass:: Forecast
	:members:
	:special-members: __init__
