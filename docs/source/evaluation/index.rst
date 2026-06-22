.. _evaluation-index:

##########
Evaluation
##########

EarthSHAB ships a small toolkit for scoring its trajectory predictions
against real flight data. The three workflows below build on each other:

* :doc:`singleEvaluation` — score one simulated flight against its APRS
  ground truth and produce a side-by-side comparison plot. Use this when
  tuning hyperparameters or debugging a single launch.
* :doc:`batchEvaluation` — run every flight in ``evaluation/launches.json``
  through the same evaluator and write a timestamped, git-tagged folder
  of results (CSV + interactive HTML summary). Use this to establish a
  baseline or measure the effect of a code change across the whole flight
  library at once.
* :doc:`batchComparison` — diff two batches (A vs B) into a single HTML
  report with per-launch deltas, win/loss counts, and overview plots.
  Use this to answer "did my change make things better, worse, or
  different?"

Every assumption (phase detection, smoothing, NaN handling, reforecast
construction, sunset detection) is shared across all three — they all
call the same :py:class:`evaluation.evaluate.BalloonEvaluator`.

.. toctree::
   :maxdepth: 1

   singleEvaluation
   batchEvaluation
   batchComparison
