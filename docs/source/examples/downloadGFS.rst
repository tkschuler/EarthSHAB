==========================
Downloading GFS Forecasts
==========================

.. important::

    As of **February 2026**, NOAA’s legacy NOMADS OpenDAP endpoint used by earlier versions of EarthSHAB has been deprecated. EarthSHAB now downloads GFS data directly using the updated NOAA distribution endpoints via GRIB access and converts them locally into NetCDF format.


Forecast data is retrieved at a **selectable spatial resolution** — ``0.25°``
(GFS ``0p25`` product, default), ``0.5°`` (``0p50``), or ``1.0°`` (``1p00``) —
and saved locally as NetCDF files (see `Unidata netcdf API <https://docs.unidata.ucar.edu/nug/current/netcdf_introduction.html>`_). This allows:

* significantly faster batch simulations
* reproducibility
* offline trajectory execution

The resolution is set with ``config`` ``_gfs_res`` (mirrored into
``netcdf_gfs['res']``). The coarser grids download far less data — a ``1.0°``
file is roughly **1/16th** the size of a ``0.25°`` file — at the cost of a more
distant nearest-grid sample (the reader picks the single nearest grid cell with
no spatial interpolation, so a coarser grid changes the resulting trajectory).

.. _gfs-filenames:

Forecast filenames
------------------

The downloaders encode both the grid resolution and the temporal step in the
output filename, so files of different grids or cadences never collide:

.. code-block:: text

   gfs_<res>_<step>_<YYYYMMDD>_<HH>.nc            e.g. gfs_0p25_3h_20260605_06.nc
   gfs_<res>_<step>_<YYYYMMDD>_<HH>_archive.nc    e.g. gfs_0p25_3h_20220822_12_archive.nc

The ``_archive`` suffix marks a file pulled from the AWS historical archive
(see :ref:`gfs-archive`). ``config`` builds this path automatically from
``_gfs_res``, ``_gfs_step_hours``, and ``forecast_start_time``, so
``forecast['file']`` and ``netcdf_gfs['nc_file']`` always point at whatever the
downloader produces.

.. note::

   Because ``forecast['file']`` is built from ``res`` and ``step_hours``, a
   stale knob can silently load a file with a different grid or cadence than
   intended. The :class:`Forecast` reader guards against this: on load it
   derives and prints the file's actual grid spacing (``resolution_deg``) and
   cadence (``resolution_hr``), then compares both against the declared
   ``netcdf_gfs['res']`` and ``netcdf_gfs['step_hours']``:

   * **GFS** files (filename is config-driven): a mismatch **raises an error** —
     it almost always means ``forecast['file']`` points at the wrong file.
   * **ERA5** files (supplied manually): a mismatch prints a non-fatal
     **warning**; the simulation proceeds on the file's actual resolution and
     cadence. ERA5 isn't produced from these GFS knobs, so the warning just
     flags that a stale config no longer describes the loaded file.


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
    * Must not be too recent (accounts for ~3–4 hour upload lag)
    * Cycles **older** than NOAA's ~9-day live window are not on NOMADS, but
      ``saveNETCDF.py`` will offer to fetch them from the AWS archive instead
      (see :ref:`gfs-archive` below).
- ``netcdf_gfs["res"]``: GFS grid resolution in degrees — ``0.25`` (default), ``0.5``, or ``1.0``. Must match the value used to build ``nc_file`` (``_gfs_res`` in ``config``); the reader enforces this (see the note under :ref:`gfs-filenames`).
- ``netcdf_gfs["lat_range"]``: Bounding-box height (in degrees) to download, centered around the ``start_coord``
- ``netcdf_gfs["lon_range"]``: Bounding-box width (in degrees) to download, centered around the ``start_coord``
- ``netcdf_gfs["download_days"]``: Determines how far into the forecast horizon to retrieve data. Maximum: ~16 days
- ``netcdf_gfs["step_hours"]``: Temporal resolution of the forecast steps to download, in hours. ``3`` (default) downloads every 3rd forecast hour; ``1`` downloads hourly. Applies to both the live and archive downloaders.

  .. important::

     Only the **0.25° grid** provides hourly steps, and only out to **f120
     (5 days)**; beyond that — and on the **0.5° and 1.0° grids at any range** —
     only 3-hourly steps exist. Requested steps that don't exist on the server
     (e.g. ``step_hours = 1`` past 120 h, or any sub-3-hourly step on a coarse
     grid) are simply skipped, so set ``step_hours = 3`` (or a multiple) for the
     0.5°/1.0° grids.

- ``netcdf_gfs["start_coord"]``: Initial latitude / longitude of the simulation. Also used as the center of the spatial subset of forecast data.

.. _gfs-archive:

Archived (historical) GFS forecasts
-----------------------------------

NOAA's live GRIB-filter service only retains roughly the last **9 days** of GFS
cycles, so historical flights cannot be re-downloaded from it. EarthSHAB also
ships ``saveNETCDF_archive.py``, which pulls the *identical* ``pgrb2.0p25``
product from the `AWS Open Data GFS archive
<https://registry.opendata.aws/noaa-gfs-bdp-pds/>`_ (bucket
``noaa-gfs-bdp-pds``, coverage from ~early 2021), using each file's ``.idx``
sidecar to fetch only the needed GRIB messages via HTTP byte-range requests. It
reads the **same** ``netcdf_gfs`` config and writes the **same v2 canonical
schema**, so an archive file is a drop-in replacement for a live download.

You do not normally call it directly. When ``saveNETCDF.py`` detects that the
requested ``forecast_start_time`` predates the live retention window, it prints a
yellow warning and asks whether to download the archived copy instead — and on
confirmation auto-switches to ``saveNETCDF_archive.py``:

.. code-block:: text

   [WARNING] Requested cycle 2022-08-22 12:00:00 is older than NOAA's ~9-day live
             retention window, so it is not available from NOMADS.
             A historical copy is available from the AWS GFS archive (noaa-gfs-bdp-pds).
             Download the archived forecast instead? [Y/n]:

In a non-interactive session (no TTY, e.g. a batch job) it prints the notice and
proceeds with the archive automatically. To run the archive downloader directly:

.. code-block:: bash

   python -m EarthSHAB.saveNETCDF_archive

Cycles **before ~2021-03** predate the AWS ``pgrb2.0p25`` archive and are
rejected; those require the `NCAR RDA ds084.1
<https://rda.ucar.edu/datasets/ds084.1/>`_ historical archive (free login).


Continuously maintaining the latest forecast (``gfs_collector.py``)
-------------------------------------------------------------------

``saveNETCDF.py`` downloads **one** configured cycle, **subset to the
``netcdf_gfs`` lat/lon box**, and exits — ideal for a small, bounded file for a
specific run. For an operational/real-time station that should always have the
freshest GFS on disk for **forward predictions**, EarthSHAB also ships
``gfs_collector.py`` — a long-running **scheduled downloader** that polls NOAA's
AWS GFS bucket and keeps a current, full-globe v2-canonical NetCDF up to date:

.. code-block:: bash

   python -m EarthSHAB.gfs_collector                    # defaults: 1.0° grid, 3-hourly steps
   python -m EarthSHAB.gfs_collector --res 0.25 --step-hours 1

On each pass it:

1. determines the latest published GFS cycle,
2. downloads any GRIB2 files for that cycle it doesn't already have (raw
   ``pgrb2`` from ``noaa-gfs-bdp-pds``, full 0–384 h horizon at the chosen step),
3. builds the v2 NetCDF **only once every server-available step for the cycle is
   on disk**, then
4. deletes the previous cycle's GRIBs and NetCDF,

and sleeps ~5 minutes before checking again.

Unlike ``saveNETCDF.py``, the collector does **not** read the ``netcdf_gfs``
config. Because it is meant to run unattended on an always-on station, it
defaults to the lightest sensible grid and cadence — **1.0° at 3-hourly steps** —
so it never silently pulls ~16× more data than intended. Override either on the
command line:

.. list-table::
   :header-rows: 1
   :widths: 20 12 68

   * - Flag
     - Default
     - Description
   * - ``--res``
     - ``1.0``
     - Grid resolution in degrees; one of ``0.25``, ``0.5``, ``1.0``.
   * - ``--step-hours``
     - ``3``
     - Hours between forecast steps. Only the ``0.25°`` grid carries
       sub-3-hourly steps, and only out to f120 (see the note above); other
       steps 404 and are skipped.

It writes the same v2 canonical schema as ``saveNETCDF.py``. Files land under
``GFS_DATA/`` at the project root (``GFS_DATA/GRIBS/`` for the raw GRIBs,
``GFS_DATA/NETCDF/`` for the output), named
``gfs_<YYYYMMDD>_<HH>z_<res>_<step>h.nc``.

.. note::

   The collector is **resilient to crashes, restarts, and partially-uploaded
   cycles**: already-downloaded GRIBs are skipped cheaply on the next pass, and
   the previous good forecast is never deleted until a *complete* new cycle has
   been converted. A truncated download therefore never replaces a good
   forecast — it just retries.

.. note::

   ``gfs_collector.py`` requires ``cfgrib`` (and its ECMWF ``eccodes`` backend)
   to read GRIB2 files. See :doc:`../installation` if it isn't already
   installed.


Running the collector in the background with ``screen``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because the collector runs forever (polling every ~5 minutes), it's best to
launch it inside a detached terminal multiplexer so it keeps running after you
log out and you can reattach to watch its output. ``screen`` is the simplest
option (``sudo apt install screen`` / ``brew install screen`` if you don't have
it). A named session called ``gfs_collector`` makes it easy to find later.

**Start it** — create the named session, activate your environment, and launch:

.. code-block:: bash

   screen -S gfs_collector              # create + attach to a session named "gfs_collector"
   conda activate earthshab             # (or your venv) — same env you run EarthSHAB in
   cd /path/to/EarthSHAB
   python -m EarthSHAB.gfs_collector    # starts polling; prints status each pass

**Detach** (leave it running in the background): press ``Ctrl-a`` then ``d``.
You'll drop back to your normal shell while the collector keeps running inside
the session.

**Check on it** — reattach to see the live output:

.. code-block:: bash

   screen -ls                           # list sessions; look for "....gfs_collector  (Detached)"
   screen -r gfs_collector              # reattach (then Ctrl-a d to detach again)

   # If it shows "(Attached)" because a stale session is holding it, force-reattach:
   screen -d -r gfs_collector

**Stop it**: reattach (``screen -r gfs_collector``) and press ``Ctrl-c`` to stop
the collector, then ``exit`` (or ``Ctrl-a`` then ``k``) to close the session.

**Respawn it** after a reboot, a crash, or an accidental close — just recreate
the session and relaunch. The collector is stateful on disk, so it picks up
right where it left off: existing GRIBs for the current cycle are skipped and
only missing/newer steps are fetched.

.. code-block:: bash

   screen -ls | grep gfs_collector || echo "not running"   # quick "is it alive?" check
   screen -S gfs_collector                                  # recreate, then activate env + relaunch as above

.. tip::

   To start it fully unattended in one line (e.g. from a startup script or over
   SSH), use a detached session and pass the command directly so you never need
   an interactive attach:

   .. code-block:: bash

      screen -dmS gfs_collector bash -lc \
        'conda activate earthshab && cd /path/to/EarthSHAB && python -m EarthSHAB.gfs_collector'

   ``tmux`` works equally well if you prefer it: ``tmux new -s gfs_collector``
   to start, ``Ctrl-b d`` to detach, ``tmux attach -t gfs_collector`` to
   reattach.