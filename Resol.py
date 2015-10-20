#!/usr/bin/python
# Time-stamp: <2015-10-20 11:02:02 marine>
## Project : Snow in the F - layer
## Subproject : functions for solving the system
## Author : Marine Lasbleis


# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import ode

#files

import param  # parameters file. Import without needed to use the param.name format
import Seismic_observations
#import eqResol
import systemdiff




def DimensionalForms(BoundConditions,PeT, Mx, Mp, Veff, gamma, X0=param.X0, d=param.d, rICp=param.rICp, rhoH=param.rhoH,alpha=param.alpha, rhoD=param.rhoD, g=param.g, hminus=param.hmin,hmaxi=param.hmax,dh=param.dt,geom="cart",AbsTole=param.AbsTol,RelTole=param.AbsTol):
    """ Integration RK (de type Dormand Prince 4-5) of the equations system for an horizontal F-layer.
    resolutionSystem(y0, t0, tmax, dt, *arg)
    Default parameters for the others arguments are :  PeT=1.e6, MX=0.17, MP=0.016, Veff=1.e7, gamma=7.2,Atol=1.e-10,Rtol=1.e-12
    The result is given in fully dimensional form.
    """
   

    if geom=="cart":
        sol = ode(systemdiff.cart).set_integrator('dopri5',atol=param.AbsTol,rtol=param.AbsTol,nsteps=1e7)
    elif geom=="spher":
        sol = ode(systemdiff.spher).set_integrator('dopri5',atol=param.AbsTol,rtol=param.AbsTol,nsteps=1e7)

    sol.set_initial_value(BoundConditions, param.hmin)
    sol.set_f_params(PeT, Mx, Mp, Veff, gamma)
    X , dXdh , phiVs ,  h = BoundConditions[0]*np.ones(1), BoundConditions[1]*np.ones(1), BoundConditions[2]*np.ones(1), param.hmin*np.ones(1)

    if param.dt>0 :
        while sol.successful() and sol.t < param.hmax:
            sol.integrate(sol.t+param.dt)
        # print("%g %g" % (r.t, r.y))
            X = np.append(X, sol.y[0])
            dXdh = np.append(dXdh, sol.y[1])
            phiVs = np.append( phiVs, sol.y[2])
            h = np.append(h, sol.t)
    else:
        while sol.successful() and sol.t > param.hmax:
            sol.integrate(sol.t+param.dt)
        # print("%g %g" % (r.t, r.y))
            X = np.append(X, sol.y[0])
            dXdh = np.append(dXdh, sol.y[1])
            phiVs = np.append( phiVs, sol.y[2])
            h = np.append(h, sol.t)
            
    
    
    #back to fully dimensional:
    x=X0*X
    z=h*d
    T=param.Tliq0*(1-Mp*h-X*Mx)

    return z,x,phiVs,T


def Vp_estimation(z,T,x,g=param.g):
    """ Estimation of the Vp profile from the results of solving the system.
    """

    DT=T-T[-1] # temperature variation in the layer compared to T[ICB]
    drhoP=-param.rhoH**2.*g*z/Seismic_observations.calcK(z)
    drhoT=-param.rhoH*param.alpha*DT#*(Mp*h+X*Mx)
    rhoL=(param.rhoD-(1-x[0])*param.rhoH-drhoT[0]-drhoP[0])/x[0]
    # print  rhoL
    ## rhoL2=x[0]/(1/(rhoD-drhoT[0]-drhoP[1])-(1-x[0])/rhoH)
    ## ## print rhoL
    
    rho_new=x*rhoL+(1-x)*param.rhoH+drhoT+drhoP
    Vp_new=np.sqrt(Seismic_observations.calcK(z)/rho_new)

    return rho_new, Vp_new
    


def SolvingEquations_DimensionalForms(BoundConditions,PeT, Mx, Mp, Veff, gamma, X0=param.X0, d=param.d, rICp=param.rICp, rhoH=param.rhoH,alpha=param.alpha, rhoD=param.rhoD, g=param.g, hminus=param.hmin,hmaxi=param.hmax,dh=param.dt,geom="cart",AbsTole=param.AbsTol,RelTole=param.AbsTol):

    # solve the equations with the given parameters
    X , dXdh , phiVs , h = eqResol.resolutionSystem(BoundConditions, hminus, hmaxi, dh, PeT, Mx, Mp, Veff, gamma, geom,Atol=AbsTole,Rtol=RelTole)

    #back to fully dimensional:
    x=X0*X
    z=h*d
    T=param.Tliq0*(1-Mp*h-X*Mx)

    return z,x,T,phiVs
