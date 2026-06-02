==========================
Downloading GFS Forecasts
==========================

.. important::

    As of **February 2026**, NOAA’s legacy NOMADS OpenDAP endpoint used by earlier versions of EarthSHAB has been deprecated. EarthSHAB now downloads GFS data directly using the updated NOAA distribution endpoints via GRIB access and converts them locally into NetCDF format.


Forecast data is retrieved at **0.25° spatial resolution** (GFS 0p25 product) and saved locally as NetCDF files (see `Unidata netcdf API <https://docs.unidata.ucar.edu/nug/current/netcdf_introduction.html>`_). This allows:

* significantly faster batch simulations
* reproducibility
* offline trajectory execution


GFS Forecast Structure
--------------------------

GFS forecasts are produced:

* **4 times per day**: 00Z, 06Z, 12Z, 18Z
* **Temporal resolution**:
  * 0–120 hours: 1-hour or 3-hour increments (depending on product)
  * 120–384 hours: 3-hour increments

Availability constraints:

* Forecasts typically become available **~3–4 hours after cycle time**
* NOAA retains approximately **9–10 days of past forecasts**
* Each forecast extends up to **384 hours (16 days)** into the future

.. warning::

   As of EarthSHAB v2.0, ``saveNETCDF.py`` writes the **v2 canonical schema**
   directly (CF-1.7, ECMWF/ERA5-style names; see :doc:`../forecast-schema-v2`).
   Forecasts downloaded with **pre-v2** EarthSHAB (the old
   `time/lev/lat/lon` + `hgtprs/tmpprs/ugrdprs/vgrdprs` layout) are no longer
   read directly — convert them once with the migration CLI (see the
   :doc:`migration guide <../migration-v2>`):

   .. code-block:: bash

      python -m EarthSHAB.forecast_processing.migrate_v1 src/EarthSHAB/forecasts/

The output NetCDF file contains:

   * Dimensions:
     * `valid_time`
     * `pressure_level` (isobaric levels, hPa, descending)
     * `latitude` (descending)
     * `longitude` (``-180`` to ``180``)

   * Variables:

     * `z` (geopotential, m²/s² — divide by g for height)
     * `t` (temperature)
     * `u` (u wind)
     * `v` (v wind)

==========================
Saving GFS Forecasts
==========================

``saveNETCDF.py`` downloads a selected GFS forecast and converts it into a v2-canonical NetCDF file for EarthSHAB. Forecasts previously downloaded from NOMADS (pre-v2) must be migrated once with ``migrate_v1`` before they can be read (see the warning above).

The forecast download is controlled via the config file.  Key parameters include:

- ``forecast_start_time``: the start time of the gfs forecast downloaded from the server. 
    * Must be one of: **00, 06, 12, 18 UTC**
    * Must be within NOAA retention window (~9 days)
    * Must not be too recent (accounts for ~3–4 hour upload lag)
- ``netcdf_gfs["lat_range"]``: How many 0.25 degree indicies to download cenetered aroudn the ``start_coord``
- ``netcdf_gfs["lon_range"]``: How many 0.25 degree indicies to download cenetered aroudn the ``start_coord``
- ``netcdf_gfs["download_days"]``: Determines how far into the forecast horizon to retrieve data. Maximum: ~16 days
- ``netcdf_gfs["start_coord"]``: Initial latitude / longitude of the simulation. Also used as the center of the spatial subset of forecast data.