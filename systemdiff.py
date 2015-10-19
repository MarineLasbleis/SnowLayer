# Time-stamp: <2015-10-16 22:56:15 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


##  Definition of the system of differential equations to be solved


def system_cart(t , y, P, Mx, Mp, V, G, Ric=6.):
    """system of equation for the spatial evolution of snow inside the F-layer
    [dy0,dy1,dy2] = system(t , y, PeT, Mx, Mp, Veff, gamma)
    """
    dy0=y[1]
    dy2=-(V+y[2])*y[1]/(y[0])
    dy1=P/Mx*(V*Mp-dy2/G*(1-Mp*t-Mx*y[0])-Mx*V*y[1])
    
    return [ dy0 , dy1 , dy2 ]


def system_spher(t , y, P, Mx, Mp, V, G, Ric=6.):
    """system of equation for the spatial evolution of snow inside the F-layer - spherical case
    [dy0,dy1,dy2] = system(t , y, PeT, Mx, Mp, Veff, gamma)
    """
    r=Ric+t
    dy0=y[1]
    dy2=-(V+y[2])*y[1]/(y[0])-2./r*y[2]
    dy1=P/Mx*(V*Mp+dy2/G*(1-Mp*t-Mx*y[0])-Mx*V*y[1])-2./r*y[1]
  
    
    return [ dy0 , dy1 , dy2 ]



def systemTest(t, y):
    return -t**2+2*t

def systemTest2(t, y, a):
    return a

def systemTest3(t , y):
    return [y[1], 4*y[0]]

def systemTest4(t , y):
    return [ 1., y[2], 4.*y[1] ]

def systemTest5(t , y, a):
    return a[0]+a[1]*t
