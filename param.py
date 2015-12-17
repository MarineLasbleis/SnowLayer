#!/usr/bin/python
# Time-stamp: <2015-12-14 17:08:23 marine>
# Project : Snow in the F - layer
# Subproject : parameters file
# Author : Marine Lasbleis

# Computation parameters

hmin = 1.
hmax = 0.
nst = 100.
dt = (hmax - hmin) / nst
AbsTol = 1.e-10
RelTol = 1.e-12
Nddy = 30


# Physical parameters

X0 = 0.1  # mass concentration
Tliq0 = 6250.  # en K
d = 200.e3  # 200km
rICp = 1221.e3
alpha = 1.95e-5  # en K-1
g = 4.4
mp = 9.e-9  # 9e-9
mX = 1.1e4  # 6000 # 1.1e4
rho0 = 12530.  # en kg.m-3
lambdaX = 1.e-9
k = 2.e-5
Vs0 = 1.e-4
rdot = 1.e-11

# Labrosse 2015
K0 = 1403.e9
Kprim0 = 3.567
Lrho = 8039.e3
Arho = 0.484
rho0 = 12451.

gam = 1.5
TL0 = 5500


rhoH = 12530
rhoD = 12060
