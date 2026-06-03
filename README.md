[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/release/python-3110/)
[![Tests](https://github.com/tkschuler/EarthSHAB/actions/workflows/tests.yml/badge.svg)](https://github.com/tkschuler/EarthSHAB/actions/workflows/tests.yml)
[![Docs](https://github.com/tkschuler/EarthSHAB/actions/workflows/docs.yml/badge.svg)](https://github.com/tkschuler/EarthSHAB/actions/workflows/docs.yml)
[![License](https://img.shields.io/badge/license-see%20LICENSE-lightgrey.svg)](LICENSE)

# EarthSHAB

> ⚠️ **v2.0.0 introduces a new unified NetCDF forecast schema.** Archived v1 forecasts
> can be converted with the `migrate_v1` CLI.
> For more information, see the **[v2 migration guide](docs/source/migration-v2.md)** (and the
> [canonical NetCDF schema](docs/source/forecast-schema-v2.md)).

Solar high altitude balloons (SHAB) are a simple and lightweight option for aerial exploration and meteorological data collection both terrestrially and on other planets. By using a
lightweight material that absorbs visual light and emits low levels of thermal radiation, solar balloons behave similarly to hot air balloons, but are capable of ascending to much higher altitudes. Unlike hot air balloons, which use a heat source to raise the temperature of the internal air, solar balloons generate heat by absorbing solar radiation, providing a free source of lift and eliminating the need for a lighter than air gas or carrying fuel.

EarthSHAB is an open source software platform for predicting the flight paths of solar balloon on Earth, adapted from [MarsSHAB](https://github.com/tkschuler/SolarBalloon), developed at the University of Arizona. Altitude profiles for a SHAB flight are generated using heat transfer modeling and dynamic analysis. By incorporating weather forecasts from NOAA, complete 3D SHAB trajectories can also be predicted.

**Features**

* 3D trajectory prediction with heat-transfer + dynamic balloon physics
* Forecast support for **GFS** (NOAA) and **ERA5** (ECMWF) reanalysis
* Multi-altitude *rainbow* trajectory predictions
* Interactive HTML trajectory and wind-vector maps
* Evaluation suite for comparing simulations against historical APRS flights, with batch / campaign reporting

See the [full documentation](https://tkschuler.github.io/EarthSHAB/) and the [Changelog](Changelog) for what's new.

## Quickstart

```
git clone https://github.com/tkschuler/EarthSHAB.git
cd EarthSHAB
pip install -r requirements.txt && pip install -e .
python -m EarthSHAB.main # run your first trajectory prediction with the included archived GFS forecast
```

See [Installation](https://tkschuler.github.io/EarthSHAB/installation.html) for the recommended conda setup (needed for `cfgrib` / ERA5).



## Examples

``config_earth.py`` includes adjustable parameters and default parameters for running any of the files discussed below. These parameters include balloon size, envelope material properties, deployment location, date and time, etc.

``saveNETCDF.py`` downloads subsets of NOAA weather forecasts for offline simulation. The temporal resolution of the download is set by ``netcdf_gfs["step_hours"]`` (3-hourly by default, ``1`` for hourly). For cycles older than NOAA's ~9-day live window it auto-switches (with a prompt) to ``saveNETCDF_archive.py``, which pulls the same forecast from the AWS GFS archive (`noaa-gfs-bdp-pds`).

``main.py*``, ``predict.py``, and ``trapezoid.py`` show examples of how to produce relevant and html-based trajectory maps using the Google maps API.


<img src = "img/rainbow_trajectories_altitude.png" />

<img src = "img/rainbow_trajectories_map.PNG" />

## Evaluation against historical flights

``evaluation/evaluate.py`` runs a single EarthSHAB simulation against a real APRS balloon track and reports how close the model came. It runs a forward simulation, a *reforecast* (forecast winds applied to the truth altitude profile to isolate wind error from vertical-motion error), detects the ascent / float / descent phases for both sim and truth, and writes a console report, a per-launch CSV, an interactive trajectory HTML, and the comparison plot below.

```
python -m evaluation.evaluate
```

<img src = "img/evaluation_comparison_SHAB14V_GFS.png" />

For sweeping a collection of flights at once and comparing across code revisions, see ``evaluation/run_batch.py`` and the [Evaluation docs](https://tkschuler.github.io/EarthSHAB/evaluation/index.html).

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for the workflow: fork the repo, branch off `main` with an issue number in the name, and open a PR. 

If EarthSHAB is useful to you, please consider **starring the repo** ⭐ to help it reach others.

## Citing EarthSHAB

If EarthSHAB played an important role in your research, then please cite the following publication
where EarthSHAB was first introduced:

[Solar Balloons - An Aerial Platform for Planetary Exploration](https://repository.arizona.edu/handle/10150/656740)
