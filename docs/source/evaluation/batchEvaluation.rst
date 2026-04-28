.. _batch-evaluation:

========================
Batch Evaluation
========================

The batch evaluation system lets you run trajectory simulations for a library of historical balloon launches and compare the results across different versions of the codebase.
Each batch is tagged with the current **git hash** so you can trace exactly which code state produced which results.

System Overview
---------------

The system lives entirely in the ``evaluation/`` directory and consists of four components:

.. list-table::
   :widths: 35 65
   :header-rows: 1

   * - File
     - Purpose
   * - ``evaluation/launches.json``
     - Master metadata file — one entry per historical launch
   * - ``evaluation/run_batch.py``
     - Runs all launches and saves results to a timestamped batch folder
   * - ``evaluation/compare_batches.py``
     - Compares metrics between two batch runs side-by-side
   * - ``tests/test_validate_launches.py``
     - Validates that all files referenced in ``launches.json`` exist and are correctly formatted


Step 1 — Populate ``launches.json``
-------------------------------------

``launches.json`` is the single source of truth for all historical flights. Each entry describes one balloon launch and points to its trajectory and forecast files.

**Required fields** (the batch will skip the launch if any are missing):

.. list-table::
   :widths: 30 15 55
   :header-rows: 1

   * - Field
     - Type
     - Description
   * - ``shab_name``
     - string
     - Balloon identifier (e.g. ``"SHAB14V"``)
   * - ``organization``
     - string
     - Launching organization (e.g. ``"UA"``)
   * - ``launch_time``
     - string
     - UTC launch time in ISO format: ``"2022-08-22 14:36:00"``
   * - ``sim_time_hr``
     - number
     - Simulation duration in hours. Set this manually — APRS trackers often lose signal near the ground at launch and landing, so auto-detection is unreliable.
   * - ``aprs_file``
     - string
     - Filename only (e.g. ``"SHAB14V-APRS.csv"``). File must exist in ``balloon_data/``.
   * - ``launch_lat``
     - number
     - Launch latitude in decimal degrees
   * - ``launch_lon``
     - number
     - Launch longitude in decimal degrees
   * - ``launch_alt_m``
     - number
     - Ground elevation at launch site in meters
   * - ``payload_weight_kg``
     - number
     - Payload mass in kg
   * - ``envelope_weight_kg``
     - number
     - Envelope mass in kg
   * - ``balloon_shape``
     - string
     - ``"sphere"`` or ``"trapezoid"``
   * - ``balloon_size``
     - number
     - Balloon diameter in meters
   * - ``gfs_file`` and/or ``era5_file``
     - string
     - Forecast filename (e.g. ``"gfs_0p25_20220822_12.nc"``). At least one is required. Set to ``null`` if not available.

**Optional fields** (fall back to current ``config_earth.py`` defaults if omitted):

``callsign``, ``landing_time``, ``areaDensityEnv``, ``cp``, ``absEnv``, ``emissEnv``, ``Upsilon``, and any ``earth_properties`` field (``Cp_air0``, ``Cv_air0``, ``Rsp_air``, ``P0``, ``emissGround``, ``albedo``).

**Example entry:**

.. code-block:: json

   {
     "shab_name": "SHAB14V",
     "organization": "UA",
     "callsign": "SHAB14V",
     "launch_time": "2022-08-22 14:36:00",
     "landing_time": null,
     "sim_time_hr": 15,
     "aprs_file": "SHAB14V-APRS.csv",
     "gfs_file": "gfs_0p25_20220822_12.nc",
     "era5_file": "SHAB14V_ERA5_20220822_20220823.nc",
     "launch_lat": 34.60,
     "launch_lon": -106.80,
     "launch_alt_m": 1000.0,
     "payload_weight_kg": 0.9,
     "envelope_weight_kg": 2.1,
     "balloon_shape": "sphere",
     "balloon_size": 5.8
   }

.. note::
   If a launch has both a GFS and an ERA5 forecast file, the batch runner will produce **two separate evaluations** — one per forecast type — in the same launch output folder.


Step 2 — Validate Before Running
----------------------------------

Before running a batch, check that all referenced files exist and are correctly formatted:

.. code-block:: bash

   pytest tests/test_validate_launches.py -v

A passing run looks like this — one row per test per launch:

.. code-block:: text

   tests/test_validate_launches.py::test_launches_json_exists PASSED
   tests/test_validate_launches.py::test_required_fields_present[UA_SHAB14V_2022-08-22] PASSED
   tests/test_validate_launches.py::test_aprs_file_exists[UA_SHAB14V_2022-08-22] PASSED
   tests/test_validate_launches.py::test_gfs_file_exists[UA_SHAB14V_2022-08-22] PASSED
   tests/test_validate_launches.py::test_era5_file_exists[UA_SHAB14V_2022-08-22] PASSED
   ...
   33 passed, 1 warning in 0.15s

The tests are **automatically parametrized** — adding a new launch to ``launches.json`` instantly adds a full set of validation tests for it.

**Extending the test suite:**

To add new validation rules, open ``tests/test_validate_launches.py`` and add a new test function. The configurable constants at the top of the file control which fields are required and what values are considered valid:

.. code-block:: python

   # To require a new field, add it here:
   REQUIRED_FIELDS = ["shab_name", "organization", ...]

   # To support a new balloon shape, add it here:
   VALID_SHAPES = {"sphere", "trapezoid"}

   # To support a new APRS CSV format, add a column set here:
   APRS_COLUMN_SETS = [
       {"time", "lat", "lng", "altitude"},             # APRS.fi standard
       {"Date", "Time(UTC)", "Latitude", "Altitude(m)"}, # Custom onboard
   ]


Step 3 — Run a Batch
----------------------

Run all launches against the **current codebase state**:

.. code-block:: bash

   python -m evaluation.run_batch --note "baseline evaluation"

The ``--note`` flag is required. Use it to describe what changed since the last batch (e.g. ``"tuned Upsilon coefficient"``). The note is stored with the results so you can remember why each batch was run.

The runner will:

1. Detect the current git hash and branch
2. Create a timestamped output folder: ``evaluation/batches/2026-04-28T1423_a3f9c12/``
3. For each launch in ``launches.json``, run a GFS and/or ERA5 evaluation
4. Skip failed launches with a full traceback printed to the console — the batch continues
5. Write a ``summary.csv`` and ``batch_info.json`` at the batch level

**Console output example:**

.. code-block:: text

   ============================================================
     EarthSHAB Batch Evaluation
     Batch ID : 2026-04-28T1423_a3f9c12
     Note     : baseline evaluation
     Launches : 2
   ============================================================

   ── UA_SHAB14V_2022-08-22 ──
     Running GFS...
     [GFS] done

     Running ERA5...
     [ERA5] done

   ── UA_SHAB1_2020-10-01 ──
     Running GFS...
     [GFS] done

   Summary → Evaluation/batches/2026-04-28T1423_a3f9c12/summary.csv
   Batch info → Evaluation/batches/2026-04-28T1423_a3f9c12/batch_info.json

   ============================================================
     Batch complete: 2/2 launches succeeded
     Total runtime: 142.3s
   ============================================================


Output Structure
-----------------

Each batch produces a self-contained folder:

.. code-block:: text

   evaluation/batches/
   └── 2026-04-28T1423_a3f9c12/
       ├── batch_info.json           ← git hash, note, runtime, launch status
       ├── summary.csv               ← all metrics, one row per launch × forecast type
       ├── UA_SHAB14V_2022-08-22/
       │   ├── SHAB14V-APRS_GFS_2022_8_22.csv
       │   ├── SHAB14V-APRS_GFS_2022_8_22.png
       │   ├── EVALUATION_SHAB14V-APRS_GFS_2022_8_22.html
       │   ├── SHAB14V-APRS_ERA5_2022_8_22.csv
       │   ├── SHAB14V-APRS_ERA5_2022_8_22.png
       │   └── EVALUATION_SHAB14V-APRS_ERA5_2022_8_22.html
       └── UA_SHAB1_2020-10-01/
           ├── SHAB1-APRS_GFS_2020_10_1.csv
           ├── SHAB1-APRS_GFS_2020_10_1.png
           └── EVALUATION_SHAB1-APRS_GFS_2020_10_1.html

**``batch_info.json``** records everything needed to reproduce or understand the batch:

.. code-block:: json

   {
     "batch_id": "2026-04-28T1423_a3f9c12",
     "note": "baseline evaluation",
     "git_hash": "a3f9c12",
     "git_branch": "devel",
     "git_commit_message": "Added batch evaluator",
     "git_dirty": false,
     "earthshab_version": "1.2.1",
     "total_runtime_s": 142.3,
     "per_launch_avg_runtime_s": 71.15,
     "launches_attempted": ["UA_SHAB14V_2022-08-22", "UA_SHAB1_2020-10-01"],
     "launches_succeeded": ["UA_SHAB14V_2022-08-22", "UA_SHAB1_2020-10-01"],
     "launches_failed": {}
   }

.. note::
   ``git_dirty: true`` means there were uncommitted changes when the batch ran. Results from dirty batches should be treated as exploratory — commit your changes before a batch you intend to keep.


Example Output
--------------

Each launch produces a **comparison plot** showing simulation vs. ground truth across three panels: altitude profile, vertical velocity (ascent / float / descent phases), and temperature/pressure.

|eval_plot|

.. |eval_plot| image:: ../../../img/evaluation_comparison_SHAB14V_GFS.png
   :width: 100%
   :alt: SHAB14V GFS evaluation comparison plot

The **metrics report** is also written as a CSV for each launch:

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


Step 4 — Compare Batches
--------------------------

After making a change to the codebase and running a second batch, compare the results:

.. code-block:: bash

   python -m evaluation.compare_batches \
       2026-04-28T1423_a3f9c12 \
       2026-04-28T1601_b7d4e21

The first argument is the **baseline** and the second is the **candidate**. The output shows metric deltas (baseline − candidate) for each launch and forecast type:

.. code-block:: text

   ========================================================================
     Batch comparison
     Baseline  (A): 2026-04-28T1423_a3f9c12
                  Note: baseline evaluation
                  Commit: Added batch evaluator
     Candidate (B): 2026-04-28T1601_b7d4e21
                  Note: tuned Upsilon coefficient
                  Commit: Tuned ascent resistance
   ========================================================================

     ── UA_SHAB14V_2022-08-22 [GFS] ──

     Metric                                         A            B      Delta (A-B)
     -------------------                   ----------   ----------   --------------
     landing distance km                      +216.697      +198.312        +18.385
     temp mae k                                +38.550       +36.210         +2.340
     sim float alt mean m                    +18988.0     +19504.0          -516.0


Workflow Summary
----------------

The recommended workflow when iterating on the codebase:

.. code-block:: bash

   # 1. Add or update launches.json entries
   # 2. Validate all data files
   pytest tests/test_validate_launches.py -v

   # 3. Run a baseline batch before making changes
   python -m evaluation.run_batch --note "before tuning Upsilon"

   # 4. Make your code changes, then run another batch
   python -m evaluation.run_batch --note "after tuning Upsilon"

   # 5. Compare
   python -m evaluation.compare_batches <batch_a_id> <batch_b_id>
