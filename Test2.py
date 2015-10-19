#!/usr/bin/python
# Time-stamp: <2015-10-19 10:19:23 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


# external libraries
import scipy.io as io
import math
import matplotlib.pyplot as plt
import numpy as np

#files
import systemdiff
import figures
import eqResol

from param import *


test_fitddy=0  ! if 1, plot figures for testing the value of ddy
figure=1
geometry="cart"
print geometry






######################
######################
####  Load seismic observations
######################
######################

Observations = io.loadmat('/Users/marine/ownCloud/Research/PREM/AK_PREM.mat')
Ks_PREM ,   alpha_AK, alpha_PREM, r_AK,  r_PREM,   rho_AK,   rho_PREM = Observations['Ks_PREM']*1e9 ,   Observations['alpha_AK']*1000., Observations['alpha_PREM'], Observations['r_AK']*1000.,  Observations['r_PREM']*1000.,   Observations['rho_AK']*1000.,  Observations['rho_PREM']

hminP=(r_PREM>=1221e3).argmax()
hmaxP=(r_PREM>=1221e3+d).argmax()
Vp_PREM=alpha_PREM[hminP+1:hmaxP+1]
radius_PREM=r_PREM[hminP+1:hmaxP+1]

hmaxA=(r_AK<=1200.5e3).argmax()
hminA=(r_AK<=1220e3+d).argmax()
Vp_AK=alpha_AK[hminA-1:hmaxA-1]
radius_AK=r_AK[hminA-1:hmaxA-1]


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

KPREM=np.array([13047, 12888, 12679, 12464])
radius=np.array([1221.5e3, 1300.0e3,1400.e3,1500.e3])


ric=np.linspace(0,3500e3,30)
Kprem_labrosse2015=(K0-K0*Kprim0*(ric**2./Lrho**2+4./5.*ric**4./Lrho**4))
rho_Labrosse2015=rho0*(1-ric**2/Lrho**2-Arho*ric**4/Lrho**4)
K_labrosse2003=1777e9*rho_Labrosse2015/7.5e3*(np.log(rho_Labrosse2015/7.5e3)+1)
K_approx=KPREM[0]-(KPREM[1]-KPREM[0])/(radius[1]-radius[0])*radius[0]+(KPREM[1]-KPREM[0])/(radius[1]-radius[0])*ric
K_approx=K_approx*1.e8

rhoH=rho0*(1-rICp**2/Lrho**2-Arho*rICp**4/Lrho**4)
drhoTad=rhoH*alpha*(-0.32e-3)*(ric-rICp)#(rhoH)**(1.-gam)*alpha*rho0**gam*TL0*2.*gam * ric**2/Lrho**2.
drhoP=-rhoH**2*g*(ric-rICp)/(KPREM[0]*1e8)
rho_test=rhoH+drhoTad+drhoP

def calcK(hicb,rICB=rICp):
    """Return the value of the bulk modulus with a linear approximation of PREM at ICB

    hicb in m
    """
    ric=rICB+hicb
    return (KPREM[0]-(KPREM[1]-KPREM[0])/(radius[1]-radius[0])*radius[0]+(KPREM[1]-KPREM[0])/(radius[1]-radius[0])*ric)*1e8


f,axa =plt.subplots(2,3)
axa[0,0].plot(r_PREM/rICp,alpha_PREM,r_AK/rICp,alpha_AK,r_PREM2/rICp,alpha_PREM2)
axa[0,0].set_title('Vp (m/s)')
axa[0,1].plot(radius_PREM/rICp,Vp_PREM,radius_AK/rICp,Vp_AK,radius_PREM2/rICp,Vp_PREM2)
axa[0,1].set_title('ZOOM - Vp (m/s)')
axa[0,2].plot(ric/rICp,Kprem_labrosse2015,r_PREM/rICp, Ks_PREM,ric/rICp,K_approx, r_PREM/rICp, alpha_PREM**2*rho_PREM, ric/rICp, K_labrosse2003, ric/rICp, 1293e9*np.ones(30))
axa[0,2].set_title('Ks (Pa)')
axa[0,2].scatter(1.,1.3047e12)
axa[0,0].axis([0,2,9000.,12000.])
axa[0,2].axis([0.5,3.,0.9e12,1.8e12])
axa[0,0].set_xlabel('r/r_icb')
axa[0,1].set_xlabel('r/r_icb')

axa[1,0].plot(r_PREM/rICp,rho_PREM,r_AK/rICp,rho_AK,ric/rICp,rho_Labrosse2015,ric/rICp,rho_test)
axa[1,0].set_title('Rho (kg/m^3)')
axa[1,1].plot(r_PREM/rICp,rho_PREM,r_AK/rICp,rho_AK,ric/rICp, rho_Labrosse2015,ric/rICp,rho_test)
axa[1,1].set_title('ZOOM - Rho (kg/m^3)')
axa[1,0].axis([0,2,11000.,13500.])
axa[1,1].axis([0.95,1.35,11800.,12200.])
axa[1,0].set_xlabel('r/r_icb')
axa[1,1].set_xlabel('r/r_icb')

print rho_test[0], rho_Labrosse2015[0]

rhoH = 12530
rhoD = 12060



######################
######################
####  Compute the density and Vp profile from solving the equations
######################
######################


def SolvingEquations_DimensionalForms(BoundConditions,PeT, Mx, Mp, Veff, gamma, X0=X0, d=d, rICp=rICp, rhoH=rhoH,alpha=alpha, rhoD=rhoD, g=g, hminus=hmin,hmaxi=hmax,dh=dt,geom="cart",AbsTole=AbsTol,RelTole=AbsTol):

    X , dXdh , phiVs , h = eqResol.resolutionSystem(BoundConditions, hminus, hmaxi, dh, PeT, Mx, Mp, Veff, gamma, geom,Atol=AbsTole,Rtol=RelTole)
    x=X0*X
    z=h*d
    T=Tliq0*(1-Mp*h-X*Mx)
    DT=T-T[-1] # temperature variation in the layer compared to T[ICB]
    drhoP=-rhoH**2.*g*z/calcK(z)
    drhoT=-rhoH*alpha*DT#*(Mp*h+X*Mx)
    rhoL=(rhoD-(1-x[0])*rhoH-drhoT[0]-drhoP[0])/x[0]
    # print  rhoL
    rhoL2=x[0]/(1/(rhoD-drhoT[0]-drhoP[1])-(1-x[0])/rhoH)
    ## print rhoL
    
    rho_new=x*rhoL+(1-x)*rhoH+drhoT+drhoP
    Vp_new=np.sqrt(calcK(z)/rho_new)
    # ax5.plot(r/1e3,Vp_new)
    
    # equation 41 Gubbins 2008 (1/rho=1/rho1+1/rho2)
    rho_new2=1./(x/rhoL2+(1-x)/rhoH)+drhoT+drhoP
    Vp_new2=np.sqrt(calcK(z)/rho_new2)
    # ax5.plot(z/1e3,Vp_new2)

    return z, x, phiVs, T, rho_new, Vp_new

###########################################
###########################################
###########################################
###########################################
###########################################
###########################################


#### Opening of the figures (if needed)

if test_fitddy :     # best value of ddy to fit the seismic values 
    fig, (ax1,ax2) = plt.subplots(1, 2, sharey=True)
    fig5, ax5 = plt.subplots()


figk, axk = plt.subplots()

#### Dimensionless parameters


for k in k*np.linspace(1.,100.,10):

    # PeT, MX, MP, Veff, gamma = 1.e6 , 0.17 , 0.016 , 1.e-7, 7.2
    PeT=Vs0*d/k
    PeX=Vs0*d/lambdaX
    Mx=X0*mX/Tliq0
    Mp=rho0*g*mp*d/Tliq0
    Veff=rdot/Vs0
    gamma=7.2

    print 'PeT {0:.3e}'.format(PeT), 'PeX {0:.3e}'.format(PeX) , 'gamma ', gamma, 'Veff {0:.3e}'.format(Veff), 'Mp {0:.3e}'.format(Mp), 'Mx {0:.3e}'.format(Mx)



    ecartVp=np.ones(Nddy) # discrpency btw Vp calculated here and Vp AK at ICB+ (to obtain the best fit)
    i=0
    DDY=np.linspace(-0.01,3,Nddy)


    for ddy in DDY:


        #### Boundary conditions at h=hmin
        BoundConditions = [1. , ddy , 0. ]

        #### Resolution of the equations
        [z, x, phiVs, T, rho_new, Vp_new]= SolvingEquations_DimensionalForms(BoundConditions,PeT, Mx, Mp, Veff, gamma)
        r=rICp+z

        #### Some plots
        if test_fitddy :
            ax1.plot(x/X0,z/1.e3)
            ax2.plot(phiVs,z/1.e3)
            ax5.plot(r/1e3,Vp_new)

        ecartVp[i]=(Vp_new[-1]-Vp_AK[-1])
        i=i+1




    a=np.argmin(ecartVp**2.)
    poly=np.polyfit(DDY[a-2:a+2],ecartVp[a-2:a+2],3)
    racines=np.roots(poly)
    ddy_OK = [racine for racine in racines if racine>DDY[a-2] and racine<DDY[a+2] ]
    ddy_OK=ddy_OK[0]




    #### Resolution of the equations with the good value for ddy
    [z, x, phiVs, T, rho_new, Vp_new]= SolvingEquations_DimensionalForms([1. , ddy_OK , 0. ],PeT, Mx, Mp, Veff, gamma, geom=geometry)

    axk.plot(z,Vp_new)

    if test_fitddy :
        ax5.plot(radius_PREM/1e3,Vp_PREM,radius_AK/1e3,Vp_AK,radius_PREM2/1e3,Vp_PREM2)
        ax5.set_title('F-layer - Vp (m/s)')
        ax5.scatter(radius_AK[-1]/1e3,Vp_AK[-1])
        #ax5.plot((radius_PREM-rICp)/1.e3,Vp_PREM,(radius_AK-rICp)/1.e3,Vp_AK)
        #ax5.scatter((radius_AK[-1]-rICp)/1.e3,Vp_AK[-1])
        ax1.set_ylabel('height above ICB (km)')
        f6,ax6=plt.subplots()
        ax6.plot(DDY,ecartVp)
        ax1.plot(x/X0,z/1.e3,'ro')
        ax2.plot(phiVs,z/1.e3,'ro')
        ax5.plot(r/1e3,Vp_new,'ro')
        ax6.scatter(ddy_OK,(Vp_new[-1]-Vp_AK[-1]))




    

    
#### Figures

if test_fitddy or figure:
    plt.show()


#figures.plotXY(X,h,'X','h')
#figures.plotXY_2(X,h, phiVs,h,0,'X','h','phiVs','h')

   
