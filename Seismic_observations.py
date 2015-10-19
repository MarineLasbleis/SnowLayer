#!/usr/bin/python
# Time-stamp: <2015-10-19 10:30:56 marine>
## Project : Snow in the F - layer
## Subproject : seismic observations
## Author : Marine Lasbleis

import math
import scipy.io as io
import numpy as np

import param



# Compute the bulk modulus K for PREM + load the PREM and AK135 profiles

def Load_PREM_AK1135_PREM2():

# Load the profiles PREM and AK135. radius in m, alpha in m/s, Ks in Pa. 
    Observations = io.loadmat('/Users/marine/ownCloud/Research/PREM/AK_PREM.mat')
    Ks_PREM ,   alpha_AK, alpha_PREM, r_AK,  r_PREM,   rho_AK,   rho_PREM = \
      Observations['Ks_PREM']*1e9 ,   Observations['alpha_AK']*1000., Observations['alpha_PREM'], Observations['r_AK']*1000., \
        Observations['r_PREM']*1000.,   Observations['rho_AK']*1000.,  Observations['rho_PREM']

    # Load the F-layer only (thickness = d, define in the param.py)
    # for PREM
    hminP=(r_PREM>=1221e3).argmax()
    hmaxP=(r_PREM>=1221e3+param.d).argmax()
    Vp_PREM=alpha_PREM[hminP+1:hmaxP+1]
    radius_PREM=r_PREM[hminP+1:hmaxP+1]
    # for AK135
    hmaxA=(r_AK<=1200.5e3).argmax()
    hminA=(r_AK<=1220e3+param.d).argmax()
    Vp_AK=alpha_AK[hminA-1:hmaxA-1]
    radius_AK=r_AK[hminA-1:hmaxA-1]

    # define the profile for PREM2 (as given in the paper)
    r_PREM2_1=np.linspace(0.,1010.0e3,10)
    alpha_PREM2_1=11.2622-6.364*(r_PREM2_1/6371.e3)**2
    r_PREM2_2=np.linspace(1010.0e3,1221.5e3,10)
    alpha_PREM2_2=11.3041-1.2730*(r_PREM2_2/6371.e3)
    r_PREM2_3=np.linspace(1221.5e3,1621.5e3,10)
    alpha_PREM2_3=4.0354+82.008*(r_PREM2_3/6371.e3)-347.769*(r_PREM2_3/6371.e3)**2+468.786*(r_PREM2_3/6371.e3)**3.
    r_PREM2=np.concatenate((r_PREM2_1,r_PREM2_2,r_PREM2_3))
    alpha_PREM2=np.concatenate((alpha_PREM2_1,alpha_PREM2_2,alpha_PREM2_3))*1000.
    radius_PREM2=r_PREM2[20:30]
    Vp_PREM2=alpha_PREM2[20:30]


    
    
