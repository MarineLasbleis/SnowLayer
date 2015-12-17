#!/usr/bin/python
# Time-stamp: <2015-12-15 18:07:04 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
import numpy as np
import json

#files

from param import *  # parameters file. Import without needed to use the param.name format
import Seismic_observations
import Resol

test_fitddy = 0  #if 1, plot figures for testing the value of ddy
figure = 1 # if 1, plot figures
geometry = "cart"  #choice of the geometry
print geometry


#Seismic_observations.Load_PREM_AK1135_PREM2()
#Seismic_observations.calcK(1.)
Seismic_observations.Figures_seism()

## Load to obtain Vp_AK(r_icb)
try:            
    with open('F_layer_AK135.dat', 'r') as file:
        data = json.load(file)
        i_ricb = data['hmax']-1
        Vp_AK = np.array(data['Vp'])[i_ricb]
        print 'Model Ak135, en r=', np.array(data['r'])[i_ricb]*1.e-3, 'km, Vp=', Vp_AK, 'km/s.'
except IOError as e:
    print "The data file does not exist. Please check for file F_layer_AK135.dat ."     

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
fig5, ax5 = plt.subplots()
figk, axk = plt.subplots()
figk2, (axk1, axk2, axk3) = plt.subplots(1, 3, sharey=True)

figdiffrho, (axdiff0, axdiff1, axdiff2) = plt.subplots(1, 3, sharey=True)

Nk = 40
K = (k*np.linspace(0.18, 1.3, Nk))# k*np.linspace(1.,1.,1)#(0.2,1.5,30)
Nmx = 40
MX = mX*np.linspace(0.15, 1.5, Nmx)
Nx = 30
XX = X0*np.linspace(0.5,1.5,Nx)

maxdiff2 = np.zeros((Nk, Nmx))

for l, mX in enumerate(MX):
# for l, X0 in enumerate(XX):
    print  l
    for i, k in enumerate(K):
    # for i, k in enumerate(K):#k*np.linspace(0.1,2.,10):
    

        PeT = Vs0*d/k
        PeX = Vs0*d/lambdaX
        Mx = X0*mX/Tliq0
        Mp = rho0*g*mp*d/Tliq0
        Veff = rdot/Vs0
        gamma = 7.2

        NDParameters = {'PeT':PeT, 'PeX':PeX, 'Mx':Mx, 'Mp':Mp, 'Veff':Veff, 'gamma':gamma} #non-dimensional parameters


        ecartVp = np.ones(Nddy)
        DDY = np.linspace(-0.01, 3, Nddy)

        for j, ddy in enumerate(DDY):  # On fait varier la condition aux limites pour ensuite trouver celle qui correspond a la valeur attendue de Vp en r=ricb. 


            BoundConditions = [1. , ddy , 0. ]

            z,x,phiVs,T = Resol.DimensionalForms(BoundConditions,NDParameters, X0, d, rICp, rhoH,
                                                 alpha, rhoD, g, hmin, hmax, dh=dt, geom="cart", AbsTole=AbsTol, RelTole=AbsTol)
            rho, Vp = Resol.Vp_estimation(z, T, x, g)
            r = rICp+z

            ax1.plot(x/X0, z/1.e3)
            ax2.plot(phiVs, z/1.e3)
            ax5.plot(r/1e3, Vp)

            ecartVp[j] = (Vp[-1]-Vp_AK)


        a = np.argmin(ecartVp**2.)
        poly = np.polyfit(DDY[a-2:a+2], ecartVp[a-2:a+2], 3)
        racines = np.roots(poly)
        ddy_OK = [racine for racine in racines if racine > DDY[a-2] and racine<DDY[a+2]]
        ddy_OK = ddy_OK[0]
        #### Resolution of the equations with the good value for ddy
        BoundConditions = [1., ddy_OK, 0. ]

        z, x, phiVs, T =  Resol.DimensionalForms(BoundConditions, NDParameters, X0, d, rICp, rhoH, \
                                                 alpha, rhoD, g, hmin, hmax, dh=dt, geom="cart", AbsTole=AbsTol, RelTole=AbsTol)
        rho, Vp = Resol.Vp_estimation(z, T, x, g)
        r = rICp+z



        ## axk.plot(z, Vp)

        ## axk1.plot(x/X0, z/1.e3, linewidth=1)
        ## axk2.plot(phiVs, z/1.e3, linewidth=1)
        ## axk3.plot(Vp, z/1.e3, linewidth=1)

        ## axdiff0.plot(rho, z/1.e3)
        ## axdiff1.plot(np.diff(rho), z[1:]/1.e3)
        ## axdiff2.plot(np.diff(rho,n=2), z[2:]/1.e3)

        maxdiff2[i,l] =  np.amax(np.diff(rho,n=2))



for fichier in ['F_layer_PREM.dat', 'F_layer_PREM2.dat', 'F_layer_AK135.dat']:
    try:
        with open(fichier, 'r') as file:
            data = json.load(file)
            hmin = data['hmin']
            hmax = data['hmax']
            radius = np.array(data['r'])[hmin:hmax]
            Vp_model = np.array(data['Vp'])[hmin:hmax]
            print file, radius[1], radius[-1]
    except IOError as e:
        print "Unable to open file", file, ". Please check the function Load_PREM_AK1135_PREM2"
    axk3.plot(Vp_model, (radius -rICp)/1.e3, linestyle='--', linewidth=2, label=fichier)
    
Seismic_observations.plotVp(z+rICp, Vp, save=1)   

axk3.set_xlim([10200, 10400])
axk3.set_ylim([0, 200])
axk1.locator_params(axis = 'x', nbins = 4)
axk2.locator_params(axis = 'x', nbins = 4)
axk3.locator_params(axis = 'x', nbins = 3)
axk1.set_title("LE concentration")
axk2.set_title("Solid flux")
axk3.set_title("Vp")
axk1.set_ylabel("Height above ICB")
axk3.legend(loc=1, prop={'size':8})

figk2.set_size_inches(8, 5)
figk2.savefig("profils.pdf", dpi=200, bbox_inches='tight', 
               transparent=True)

fi, a = plt.subplots()
a.plot(K*10400.*750., -maxdiff2)
fi.savefig("diff2.pdf", dpi=200, bbox_inches='tight', 
               transparent=True)

X, Y = np.meshgrid(MX/1e3, K*10400*750)
# X, Y = np.meshgrid(XX, K*10400*750)
fi2, a2 = plt.subplots()
mat = a2.contourf(X, Y, -maxdiff2, 30)
mat2 = a2.contour(X, Y, -maxdiff2>0, 1)
fi2.colorbar(mat)
a2.set_xlabel('mX (1e3 K)')
a2.set_ylabel('k (W/K/m)')
fi2.savefig("diff2_matrix.pdf", dpi=200, bbox_inches='tight', 
               transparent=True)

# plt.show()

