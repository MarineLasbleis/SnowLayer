#!/usr/bin/python
# Time-stamp: <2015-12-13 20:06:19 marine>
# Project : Snow in the F - layer
# Subproject : seismic observations
# Author : Marine Lasbleis

import math
import scipy.io as io
import numpy as np
import json  # for writing in files
import os.path
import matplotlib.pyplot as plt

import param


# Compute the bulk modulus K for PREM + load the PREM and AK135 profiles

def Load_PREM_AK135_PREM2():

    # Load the profiles PREM and AK135. radius in m, alpha in m/s, Ks in Pa.
    Observations = io.loadmat(
        '/Users/marine/ownCloud/Research/PREM/AK_PREM.mat')
    Ks_PREM, alpha_AK, alpha_PREM, r_AK, r_PREM, rho_AK, rho_PREM = \
        Observations['Ks_PREM'] * 1e9,   Observations['alpha_AK'] * 1000., Observations['alpha_PREM'], Observations['r_AK'] * 1000., \
        Observations['r_PREM'] * 1000.,   Observations['rho_AK'] * \
        1000.,  Observations['rho_PREM']

    # Load the F-layer only (thickness = d, define in the param.py)
    # for PREM
    hminP = (r_PREM >= 1221e3).argmax()
    hmaxP = (r_PREM >= 1221e3 + param.d).argmax()
    Vp_PREM = alpha_PREM[hminP + 1:hmaxP + 1]
    radius_PREM = r_PREM[hminP + 1:hmaxP + 1]
    # for AK135
    hmaxA = (r_AK <= 1200.5e3).argmax()
    hminA = (r_AK <= 1220e3 + param.d).argmax()
    # Vp_AK=alpha_AK[hminA-1:hmaxA-1]
    # radius_AK=r_AK[hminA-1:hmaxA-1]
    K_AK = np.zeros(r_AK.size)

    # define the profile for PREM2 (as given in the paper)
    r_PREM2_1 = np.linspace(0., 1010.0e3, 10)
    alpha_PREM2_1 = 11.2622 - 6.364 * (r_PREM2_1 / 6371.e3)**2
    r_PREM2_2 = np.linspace(1010.0e3, 1221.5e3, 10)
    alpha_PREM2_2 = 11.3041 - 1.2730 * (r_PREM2_2 / 6371.e3)
    r_PREM2_3 = np.linspace(1221.5e3, 1621.5e3, 10)
    alpha_PREM2_3 = 4.0354 + 82.008 * (r_PREM2_3 / 6371.e3) - 347.769 * (
        r_PREM2_3 / 6371.e3)**2 + 468.786 * (r_PREM2_3 / 6371.e3)**3.
    r_PREM2 = np.concatenate((r_PREM2_1, r_PREM2_2, r_PREM2_3))
    alpha_PREM2 = np.concatenate(
        (alpha_PREM2_1, alpha_PREM2_2, alpha_PREM2_3)) * 1000.
    # radius_PREM2=r_PREM2[20:30]
    # Vp_PREM2=alpha_PREM2[20:30]
    K_PREM2 = np.zeros(30)
    rho_PREM2 = np.zeros(30)

    # Bulk modulus
    # Bulk modulus K: from PREM
    KPREM = Ks_PREM
    # from Labrosse 2003/2015
    ric = np.linspace(0, 3500e3, 30)
    Kprem_labrosse2015 = (param.K0 - param.K0 * param.Kprim0 *
                          (ric**2. / param.Lrho**2 + 4. / 5. * ric**4. / param.Lrho**4))
    rho_Labrosse2015 = param.rho0 * \
        (1 - ric**2 / param.Lrho**2 - param.Arho * ric**4 / param.Lrho**4)
    K_labrosse2003 = 1777e9 * rho_Labrosse2015 / \
        7.5e3 * (np.log(rho_Labrosse2015 / 7.5e3) + 1)
    K_approx = KPREM[0] - (KPREM[1] - KPREM[0]) / (radius_PREM[1] - radius_PREM[0]) * \
        radius_PREM[0] + (KPREM[1] - KPREM[0]) / \
        (radius_PREM[1] - radius_PREM[0]) * ric
    K_approx = K_approx * 1.e8

    DATA_PREM = {'r': r_PREM.tolist(), 'Vp': alpha_PREM.tolist(), 'K': KPREM.tolist(
    ), 'hmin': hminP + 1, 'hmax': hmaxP + 1, 'rho': rho_PREM.tolist()}
    f = open("F_layer_PREM.dat", 'w')
    json.dump(DATA_PREM, f)
    f.close()
    DATA_PREM2 = {'r': r_PREM2.tolist(), 'Vp': alpha_PREM2.tolist(
    ), 'K': K_PREM2.tolist(), 'hmin': 20, 'hmax': 30, 'rho': rho_PREM2.tolist()}
    f = open("F_layer_PREM2.dat", 'w')
    json.dump(DATA_PREM2, f)
    f.close()
    DATA_AK135 = {'r': r_AK.tolist(), 'Vp': alpha_AK.tolist(
    ), 'hmin': hminA - 1, 'K': K_AK.tolist(), 'hmax': hmaxA - 1, 'rho': rho_AK.tolist()}
    f = open("F_layer_AK135.dat", 'w')
    json.dump(DATA_AK135, f)
    f.close()
    Vp_Labrosse = np.sqrt(Kprem_labrosse2015 / rho_Labrosse2015)
    DATA_Labrosse = {'r': ric.tolist(), 'Vp': Vp_Labrosse.tolist(
    ), 'hmin': 1, 'K': Kprem_labrosse2015.tolist(), 'hmax': 30, 'rho': rho_Labrosse2015.tolist()}
    f = open("F_layer_Labrosse.dat", 'w')
    json.dump(DATA_Labrosse, f)
    f.close()


def calcK(hicb, rICB=param.rICp):
    """Return the value of the bulk modulus with a linear approximation of PREM at ICB

    hicb in m, K in Pa.
    
    """
    if os.path.isfile('F_layer_PREM.dat') == 0:
        print "Load seismic informations."
        print ".dat files containing seismic informations are created."
        Load_PREM_AK135_PREM2()
    try:
        with open('F_layer_PREM.dat', 'r') as file:
            data = json.load(file)
            hmin = data['hmin']
            hmax = data['hmax']
            radius = np.array(data['r'])
            radius = radius[hmin:hmax]
            K = np.array(data['K'])
            K = K[hmin:hmax]
    except IOError as e:
        print "Unable to open file DATA_PREM.dat. Please check the function Load_PREM_AK1135_PREM2"

    ric = rICB + hicb
    return (K[0] - (K[1] - K[0]) / (radius[1] - radius[0]) * radius[0] + (K[1] - K[0]) / (radius[1] - radius[0]) * ric)


def Figures_seism():
    """ plot seismic observations for the F-layer: PREM, PREM2, AK135 and fit by Labrosse 2015.

    Use files : 'F_layer_PREM.dat', 'F_layer_PREM2.dat', 'F_layer_AK135.dat', 'F_layer_Labrosse.dat'
    Data are json files with: r, Vp, K, rho, hmax, hmin
    If a file does not exist, please run Load_PREM_AK135_PREM2() to obtain data from the Matlab file.
    6 figures :
    0 1 2
    3 4 5
    0: Vp in the whole Earth
    1: Vp zoom in the F-layer
    2: K in the F-layer
    3: density profile
    
    """
    
    f, axa = plt.subplots(2, 3)
    #  0 1 2
    #  3 4 5
    # 0: Vp in the whole Earth
    # 1: Vp zoom in the F-layer
    # 2: K in the F-layer
    # 3: density profile

    for fichier in ['F_layer_PREM.dat', 'F_layer_PREM2.dat', 'F_layer_AK135.dat', 'F_layer_Labrosse.dat']:

        try:
            with open(fichier, 'r') as file:
                print "=========", fichier, "========="
                data = json.load(file)
                r = np.array(data['r']) / param.rICp
                Vp = np.array(data['Vp'])
                K = np.array(data['K'])
                rho = np.array(data['rho'])
                hmax = data['hmax']
                hmin = data['hmin']

                axa[0, 0].plot(r, Vp, label=fichier)
                axa[0, 1].plot(r[hmin:hmax], Vp[hmin:hmax], label=fichier)
                axa[0, 2].plot(r, K, label=fichier)
                axa[1, 0].plot(r, rho, label=fichier)
                axa[1, 1].plot(r, rho, label=fichier)
               # axa[1,0].plot(r,rho)

        except IOError as e:
            print "The data file does not exist. Please check the function Load_PREM_AK1135_PREM2 and particularly for file ", files

    axa[0, 0].set_title('Vp (m/s)')
    axa[0, 0].axis([0, 2, 9000., 12000.])
    axa[0, 0].set_xlabel('r/r_icb')
    axa[0, 1].set_title('ZOOM - Vp (m/s)')
    axa[0, 1].set_xlabel('r/r_icb')
    axa[0, 1].axis([0.95, 1.35, 10200, 10400])
    axa[0, 2].set_title('Ks (Pa)')
    axa[0, 2].scatter(1., 1.3047e12)
    axa[0, 2].axis([0., 2., 0.9e12, 1.8e12])
    axa[1, 0].set_title('Rho (kg/m^3)')
    axa[1, 0].axis([0, 2, 11000., 13500.])
    axa[1, 0].set_xlabel('r/r_icb')
    axa[1, 1].set_title('ZOOM - Rho (kg/m^3)')
    axa[1, 1].axis([0.95, 1.35, 11800., 12200.])
    axa[1, 1].set_xlabel('r/r_icb')
    axa[0, 0].legend(prop={'size': 10})
    axa[0, 1].legend(prop={'size': 10})
    axa[0, 2].legend(prop={'size': 8})
    axa[1, 0].legend(prop={'size': 8})
    axa[1, 1].legend(prop={'size': 8})

    # plt.show()


def plotVp(z_, Vp_, save=0, name="Vp.pdf"):

    fig, ax = plt.subplots()
    for fichier in ['F_layer_PREM.dat', 'F_layer_PREM2.dat', 'F_layer_AK135.dat']:
        try:
            with open(fichier, 'r') as file:
                data = json.load(file)
                hmin = data['hmin']
                hmax = data['hmax']
                radius = np.array(data['r'])[hmin:hmax]
                Vp = np.array(data['Vp'])[hmin:hmax]
                print file, radius[1], radius[-1]
        except IOError as e:
            print "Unable to open file", file, ". Please check the function Load_PREM_AK1135_PREM2"
        ax.plot(Vp / 1e3, radius / param.rICp)
    ax.plot(Vp_ / 1e3, z_ / param.rICp)
    if save == 1:
        plt.savefig(name)

if __name__ == '__main__':

    Load_PREM_AK135_PREM2()
    Figures_seism()

    plt.show()
