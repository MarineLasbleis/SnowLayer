#!/usr/bin/python
# Time-stamp: <2015-11-05 22:31:11 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
import numpy as np
import json

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
Seismic_observations.Figures_seism()

## Load to obtain Vp_AK(r_icb)
try:            
    with open('F_layer_AK135.dat','r') as file:
        data=json.load(file)
        i_ricb=data['hmax']-1
        Vp_AK=np.array(data['Vp'])[i_ricb]
        print 'Model Ak135, en r=', np.array(data['r'])[i_ricb]*1.e-3, 'km, Vp=', Vp_AK, 'km/s.'
except IOError as e:
    print "The data file does not exist. Please check for file F_layer_AK135.dat ."     

fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)
fig5, ax5 = plt.subplots()
figk, axk = plt.subplots()
figk2, (axk1,axk2) = plt.subplots(1, 2, sharey=True)


for k in k*np.linspace(0.1,2.,10):

    PeT=Vs0*d/k
    PeX=Vs0*d/lambdaX
    Mx=X0*mX/Tliq0
    Mp=rho0*g*mp*d/Tliq0
    Veff=rdot/Vs0
    gamma=7.2

    NDParameters={'PeT':PeT, 'PeX':PeX, 'Mx':Mx, 'Mp':Mp, 'Veff':Veff, 'gamma':gamma} #non-dimensional parameters

    i=0
    ecartVp=np.ones(Nddy)
    DDY=np.linspace(-0.01,3,Nddy)

    for ddy in DDY:  # On fait varier la condition aux limites pour ensuite trouver celle qui correspond a la valeur attendue de Vp en r=ricb. 


        BoundConditions = [1. , ddy , 0. ]

        z,x,phiVs,T = Resol.DimensionalForms(BoundConditions,NDParameters, X0, d, rICp, rhoH,alpha, rhoD, g, hmin,hmax,dh=dt,geom="cart",AbsTole=AbsTol,RelTole=AbsTol)
        rho, Vp = Resol.Vp_estimation(z,T,x,g)
        r=rICp+z

        ax1.plot(x/X0,z/1.e3)
        ax2.plot(phiVs,z/1.e3)
        ax5.plot(r/1e3,Vp)

        ecartVp[i], i=(Vp[-1]-Vp_AK), i+1


    a=np.argmin(ecartVp**2.)
    poly=np.polyfit(DDY[a-2:a+2],ecartVp[a-2:a+2],3)
    racines=np.roots(poly)
    ddy_OK = [racine for racine in racines if racine>DDY[a-2] and racine<DDY[a+2] ]
    ddy_OK=ddy_OK[0]
    #### Resolution of the equations with the good value for ddy
    BoundConditions = [1. , ddy_OK , 0. ]

    z,x,phiVs,T=  Resol.DimensionalForms(BoundConditions,NDParameters, X0, d, rICp, rhoH,alpha, rhoD, g, hmin,hmax,dh=dt,geom="cart",AbsTole=AbsTol,RelTole=AbsTol)
    rho, Vp = Resol.Vp_estimation(z,T,x,g)
    r=rICp+z
    axk.plot(z,Vp)

    axk1.plot(x/X0,z/1.e3)
    axk2.plot(phiVs,z/1.e3)

Seismic_observations.plotVp(z+rICp,Vp)   



plt.show()
