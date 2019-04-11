# Daypym
Daypym is a project to preprocess and run DAYSIM solar irradiance simulations with Python. It uses EnergyPlus IDF files, as geometry input.

# :warning:Deprecated!
**You can find the latest version of Daypym here:**
[https://gitlab.tue.nl/bp-tue/daypym](https://gitlab.tue.nl/bp-tue/daypym) :heavy_check_mark:

# Features

## Geometry preprocessing

* Translating [EnergyPlus](https://energyplus.net/) IDF geometry to .rad file for [DAYSIM](http://daysim.ning.com/) simulation (utilizing [Eppy](https://github.com/santoshphilip/eppy) and [GeomEppy](https://github.com/jamiebull1/geomeppy))
* Generating sensor points over selected surfaces (utilizing GeomEppy)

## Running DAYSIM simulations

* Creating folder structure for the simulations
* Writing .hea file for DAYSIM
* Running DAYSIM simulations

## Postprocessing

* Visualizing irradiance simulation results

# "Installation"

* Install DAYSIM
* Install Eppy
* Install GeomEppy
* Then:
```python
import sys
sys.path.append(r'Path/Of/Folder/Containing/Daypym.py')
from daypym import *
```
