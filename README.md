[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)

# EarthSHAB

Solar high altitude balloons (SHAB) are a simple and lightweight option for aerial exploration and meteorological data collection both terrestrially and on other planets. By using a
lightweight material that absorbs visual light and emits low levels of thermal radiation, solar balloons behave similarly to hot air balloons, but are capable of ascending to much higher altitudes. Unlike hot air balloons, which use a heat source to raise the temperature of the internal air, solar balloons generate heat by absorbing solar radiation, providing a free source of lift and eliminating the need for a lighter than air gas or carrying fuel.

EarthSHAB is an open source software platform for predicting the flight paths of solar balloon on Earth, adapted from [MarsSHAB](https://github.com/tkschuler/SolarBalloon), developed at the University of Arizona. Altitude profiles for a SHAB flight are generated using heat transfer modeling and dynamic analysis. By incorporating weather forecasts from NOAA, complete 3D SHAB trajectories can also be predicted.  

## Installation

This simulation has been tested to run on Ubuntu 22.04 and Python 3.9.


For easy install on WSL and Ubuntu use:
```
pip install -r requirements.txt
```

Setup Environment:
```
pip install -e .
```

## Examples

``config_earth.py`` includes adjustable parameters and default parameters for running any of the files discussed below. These parameters include balloon size, envelope material properties, deployment location, date and time, etc.

``saveNETCDF.py`` downloads subsets of NOAA weather forecasts for offline simulation

``main.py*``, ``predict.py``, and ``trapezoid.py`` show examples of how to produce relevant and html-based trajectory maps using the Google maps API.

These examples can all be run with the included GFS and ERA forecasts as well as a SHAB balloon trajectory (SHAB14-V) in the required APRS.fi csv format.

<img src = "img/rainbow_trajectories_altitude.png" />

<img src = "img/rainbow_trajectories_map.PNG" />


## Author

* **Tristan Schuler** - *U.S. Naval Research Laboratory, University of Arizona*
* Additional Contributions: **Craig Motell** - *NIWC Pacific*

## Acknowledgments

Hat Tip to [Bovine Aerospace](https://bovineaerospace.wordpress.com/), who developed an initial solar balloon model in C++. This code was adapted from their [repo](https://github.com/tunawhiskers/balloon_trajectory).

## Citing EarthSHAB

If EarthSHAB played an important role in your research, then please cite the following publication
where EarthSHAB was first introduced:

[Solar Balloons - An Aerial Platform for Planetary Exploration](https://repository.arizona.edu/handle/10150/656740)
