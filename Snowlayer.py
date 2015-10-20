#!/usr/bin/python
# Time-stamp: <2015-10-20 11:03:58 marine>
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
import Resol

test_fitddy=0  #if 1, plot figures for testing the value of ddy
figure=1 # if 1, plot figures
geometry="cart"  #choice of the geometry
print geometry


#Seismic_observations.Load_PREM_AK1135_PREM2()
#Seismic_observations.calcK(1.)
#Seismic_observations.Figures_seism()



PeT=Vs0*d/k
PeX=Vs0*d/lambdaX
Mx=X0*mX/Tliq0
Mp=rho0*g*mp*d/Tliq0
Veff=rdot/Vs0
gamma=7.2
BoundConditions = [1. , 0.4 , 0. ]

z,x,phiVs,T = Resol.DimensionalForms(BoundConditions,PeT, Mx, Mp, Veff, gamma, X0, d, rICp, rhoH,alpha, rhoD, g, hmin,hmax,dh=dt,geom="cart",AbsTole=AbsTol,RelTole=AbsTol)
rho, Vp = Resol.Vp_estimation(z,T,x,g)

r=rICp+z

fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)
fig5, ax5 = plt.subplots()
ax1.plot(x/X0,z/1.e3)
ax2.plot(phiVs,z/1.e3)
ax5.plot(r/1e3,Vp)

plt.show()
