#!/usr/bin/python
# Time-stamp: <2015-10-19 10:26:22 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
import numpy as np

#files

from param import *  # parameters file. Import without needed to use the param.name format
import Seismic_observations


test_fitddy=0  #if 1, plot figures for testing the value of ddy
figure=1 # if 1, plot figures
geometry="cart"  #choice of the geometry
print geometry


Seismic_observations.Load_PREM_AK1135_PREM2()
