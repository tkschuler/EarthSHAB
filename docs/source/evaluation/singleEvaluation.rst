.. _single-evaluation:

====================
Single Evaluation
====================

A single evaluation compares one simulation run against a historical balloon flight.
It uses the current ``config_earth.py`` directly — no metadata file or batch runner needed.

Prerequisites
-------------

You need two things before running an evaluation:

1. **A downloaded forecast file** (GFS or ERA5) in ``src/EarthSHAB/forecasts/``
2. **A balloon trajectory CSV** in ``src/EarthSHAB/balloon_data/`` in the `APRS.fi <https://aprs.fi>`_ format

.. note::
   Balloon trajectory data can be downloaded after landing from `APRS.fi <https://aprs.fi>`_.
   Search for your callsign and export the track as CSV.


Configure ``config_earth.py``
------------------------------

Open ``src/EarthSHAB/config_earth.py`` and update the following sections to match your flight:

**1. Point to your APRS trajectory file:**

.. code-block:: python

   balloon_trajectory = parent_dir + "balloon_data/SHAB14V-APRS.csv"

**2. Set the simulation start time to match launch:**

.. code-block:: python

   start_time = datetime.fromisoformat("2022-08-22 14:36:00")  # UTC

**3. Set the launch coordinates:**

.. code-block:: python

   simulation = dict(
       start_time = start_time,
       sim_time = 15,        # hours — set 1–2 hours beyond actual flight duration
       start_coord = {
           "lat": 34.60,     # launch latitude
           "lon": -106.80,   # launch longitude
           "alt": 1000.,     # ground elevation (m)
           "timestamp": start_time,
       },
       min_alt = 1000.,
       balloon_trajectory = balloon_trajectory,
       ...
   )

**4. Set balloon physical properties:**

.. code-block:: python

   balloon_properties = dict(
       shape = 'sphere',
       d = 5.8,        # diameter (m)
       mp = 0.9,       # payload mass (kg)
       mEnv = 2.1,     # envelope mass (kg)
       ...
   )

**5. Select your forecast type and file:**

.. code-block:: python

   forecast = dict(
       forecast_type = "GFS",   # or "ERA5"
       forecast_start_time = "2022-08-22 12:00:00",
       GFSrate = 60,
   )


Run the Evaluation
------------------

From the repository root:

.. code-block:: bash

   python -m evaluation.evaluate

The evaluator will:

1. Run the **forward simulation** using the forecast winds
2. Run the **reforecast simulation** — forecast winds applied to the actual APRS altitude profile
3. Compare both against the ground-truth APRS track
4. Print a metrics report to the console
5. Save outputs to ``evaluation/``


Single Evaluation Output
------------------------

**Console report:**

.. code-block:: text

   ====================================================================
     EarthSHAB Evaluation: SHAB14V-APRS
   ====================================================================
     Metric                               Sim      Truth       Diff  Unit
   --------------------------------------------------------------------
     Float Alt Mean (m)                 18988      20396      -1407  m
     Float Alt Std (m)                    378        447        -69  m
     Ascent Rate Mean (m/s)              1.91       2.36      -0.45  m/s
     Ascent Rate Std (m/s)               0.33       0.96      -0.63  m/s
     Descent Rate Mean (m/s)            -2.82      -2.37      -0.44  m/s
     Elapsed Time (min)                 779.8      789.0       -9.3  min
     Landing Lat (°)                  34.7816    34.5468     0.2347  °
     Landing Lon (°)                -106.7869  -109.1340     2.3471  °
   --------------------------------------------------------------------
     Distance Off (m)                                        216697  m
     Landing Time (MST)  2022-08-22 20:35   2022-08-22 20:46  -11.2  min
   --------------------------------------------------------------------
     Temperature MAE                                          38.55  K
     Pressure MAE                                               659  Pa
   --------------------------------------------------------------------
     GFS Forecast + Truth Altitude (reforecast landing vs truth)
     Distance Off (m)                                         44811  m
   ====================================================================

   Start-time analysis
   Current config start_time : 2022-08-22 14:36:00 UTC
   First APRS transmission   : 2022-08-22 14:37:53 UTC  (1557 m)
   Estimated launch time     : 2022-08-22 14:32:23 UTC
   Suggested start_time      : "2022-08-22 14:32:23"

.. tip::
   The **start-time analysis** at the bottom extrapolates the initial APRS ascent rate back to the
   launch elevation and suggests a more accurate ``start_time``. If the suggested time differs
   significantly from your configured value, update ``config_earth.py`` and re-run.

**Comparison plot** — three panels showing altitude profile, vertical velocity by phase, and temperature/pressure vs. ISA model:

|eval_single_plot|

.. |eval_single_plot| image:: ../../../img/evaluation_comparison_SHAB14V_GFS.png
   :width: 100%
   :alt: SHAB14V single evaluation comparison plot (GFS)

**Saved files** (written to ``evaluation/``):

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - File
     - Contents
   * - ``SHAB14V-APRS_GFS_2022_8_22.csv``
     - All metrics in tabular form (Sim, Truth, Diff)
   * - ``SHAB14V-APRS_GFS_2022_8_22.png``
     - Comparison plot (three-panel figure)
   * - ``EVALUATION_SHAB14V-APRS_GFS_2022_8_22.html``
     - Interactive trajectory map


Interpreting Results
--------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Metric
     - What it tells you
   * - **Float Alt Mean**
     - How close the predicted float altitude is to actual. A large negative diff means the sim floated lower than reality — check envelope or payload mass.
   * - **Ascent Rate Mean**
     - Differences here often indicate incorrect balloon diameter or envelope properties.
   * - **Distance Off**
     - Great-circle distance between simulated and actual landing. Primarily driven by forecast wind accuracy.
   * - **Reforecast Distance Off**
     - Landing error when the *actual* altitude profile is used with forecast winds. Isolates horizontal wind error from vertical motion error.
   * - **Temperature MAE**
     - Difference between simulated atmospheric temperature and onboard sensor. Large values suggest radiation model tuning may help.
