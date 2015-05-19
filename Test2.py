# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
import numpy as np

#files
import systemdiff
import figures
import eqResol





#### Computation parameters

hmin=1.
hmax=0.
nst=1000.;
dt=(hmax-hmin)/nst
AbsTol=1.e-10
RelTol=1.e-12

#### Physical parameters

X0=0.08  # mass concentration
Tliq0=6250. # en K
d=200000. # 200km
alpha=1.95e-5 #en K-1
g=4.4
mp=9.e-9 # 9e-9
mX=1.e4 # 1.1e4
rho0=12530. #en kg.m-3
lambdaX=1.e-9
k=1.8e-5
Vs0=1.e-4
rdot=1.e-11

#### Dimensionless parameters

# PeT, MX, MP, Veff, gamma = 1.e6 , 0.17 , 0.016 , 1.e-7, 7.2
PeT=Vs0*d/k
PeX=Vs0*d/lambdaX
Mx=X0*mX/Tliq0
Mp=rho0*g*mp*d/Tliq0
Veff=rdot/Vs0
gamma=7.2




print 'PeT {0:.3e}'.format(PeT), 'PeX {0:.3e}'.format(PeX) , 'gamma ', gamma, 'Veff {0:.3e}'.format(Veff), 'Mp {0:.3e}'.format(Mp), 'Mx {0:.3e}'.format(Mx)


####  Load seismic observations
Observations = io.loadmat('/Users/marine/ownCloud/Research/PREM/PREM+AK_pourVp.mat')
Ks_PREM ,   alpha_AK, alpha_PREM, r_AK,  r_PREM,   rho_AK,   rho_PREM = Observations['Ks_PREM'] ,   Observations['alpha_AK'], Observations['alpha_PREM'], Observations['r_AK']*1000,  Observations['r_PREM']*1000,   Observations['rho_AK'],  Observations['rho_PREM']

hminP=(r_PREM>=1221e3).argmax()
hmaxP=(r_PREM>=1221e3+d).argmax()
Vp_PREM=alpha_PREM[hminP:hmaxP]
radius_PREM=r_PREM[hminP:hmaxP]

hmaxA=(r_AK<=1217.5e3).argmax()
hminA=(r_AK<=1220e3+d).argmax()
Vp_AK=alpha_AK[hminA:hmaxA]
radius_AK=r_AK[hminP:hmaxP]




fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)

for ddy in np.linspace(-0.5,2,20):

    #### Boundary conditions at h=hmin
    BoundConditions = [1. , ddy , 0. ]



    #### Resolution of the equations

    X , dXdt , phiVs , h = eqResol.resolutionSystem3(BoundConditions, hmin, hmax, dt, PeT, Mx, Mp, Veff, gamma,Atol=AbsTol,Rtol=RelTol)

   
    ax1.plot(X,h)
    ax2.plot(phiVs,h)


plt.show()
    
#### Figures

#figures.plotXY(X,h,'X','h')
#figures.plotXY_2(X,h, phiVs,h,0,'X','h','phiVs','h')

   
