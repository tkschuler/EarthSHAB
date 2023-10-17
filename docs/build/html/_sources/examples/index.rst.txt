.. _examples-index:

###################
EarthSHAB Usage
###################

.. toctree::
   :maxdepth: 1
   :caption: Downloading and Saving Forecasts:
   
   downloadGFS
   downloadERA5

These files give examples on how to use EarthSHAB. 

Configuration File
===================
 

``config.py`` is the heart of any EarthSHAB simulation and the only file that needs to be updated to use EarthSHAB at the current release.  Many parameters can be adjusted for running EarthSHAB simulations and the parameters are grouped together into convienent categories: 

- **balloon_properties:** includes parameters for adjusting the balloon size and material propeties.  Currently only the 'sphere' shape is supported.  If using a standard 6m charcoal-coated SHAB balloon, the payload and envelope masses are the only parameters that need to be updated. 
- **netcdf_gfs:** includes parameters for reading and saving GFS netcdf forecasts from `NOAA's NOMADS server <https://nomads.ncep.noaa.gov/>`_. Currently only supporting forecasts with a resolution of 0.25 degrees for better trajectroy prediction accuracy. 
- **netcdf_era5:** 
- **simulation:** 
- **forecast_type:** Either GFS or ERA5.  For either to work, a corresponding netcdf file must be download and in the forecasts directory. 
- **GFS.GFSRate:** adjusts how often wind speeds are looked up (default is once per 60 seconds). The fastest look up rate is up to 1.0 s, which provides the highest fidelity predictions, but also takes significantly longer to run a simulation.  
- **dt:** is default to 1.0s integrating intervals.
- **earth_properties:** include mostly constant proeprties for earth.  These values typically don't need to be changed, however the ground emissivity and albedo might be changed if predicting balloon launches over oceans or snow covered locations. 


.. tip:: If when increasing ``dt``, nans are introduced,  this is a numerical integration issue that can be solved by decreasing dt. Most values over 2.5s cause this issue when running a simulation. 



Examples
===================

``main.py`` performs a full-physics based simulation taking into account: balloon geometry and properties, radiation sources for the coordinate and altitude, dynamic force analysis, and weather forecasts.  The simulation outputs an altitude plot, balloon temperature, a 3d wind rose and atmospheric temperature profile for the starting location, and an interactive google maps trajectory. 

**New** in v1.1, ``main.py`` supports comparing simulations (using either a GFS forecast or ERA5 renalysis) to a historical balloon flight. At this time, only balloon trajectories in the *aprs.fi* format are supported. Flights using APRS devices for real time tracking can be downloaded after landing from  `APRS.fi <https:aprs.fi>`_.  To use this functionality, change ``balloon_trajectory`` in the configuration file from None to a filename.

|pic0|

.. |pic0| image:: ../../../img/SHAB12-V.png
   :width: 100%


``predict.py`` uses EarthSHAB to generate a family of altitude and trajectory predictions at different float altitudes. Estimating float altitudes for a standard charcoal SHAB balloon are challenging and can range from SHAB to SHAB even if all the **balloon_properties** are the same and the balloons are launches on the same day. This example produces a rainbow altitude plot and interactive google maps trajectory prediction. The float altitudes are adjusted by changing the payload mass in .25 increments to simulate varrying float altitudes.

|pic1| |pic2|

.. |pic1| image:: ../../../img/rainbow_trajectories_altitude.png
   :width: 48%

.. |pic2| image:: ../../../img/rainbow_trajectories_map.PNG
   :width: 48%

``trajectory.py`` provides a rough example of how to use EarthSHAB with a manual altitude profile to estimate wind-based trajectory predictions.  Alternatively, you can generate a csv file in the *aprs.fi* format and inclue it in ``config.py``. 


