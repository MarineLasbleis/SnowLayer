## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


##  Definition of the system of differential equations to be solved


def system(t , y, P, Mx, Mp, V, G):
    """system of equation for the spatial evolution of snow inside the F-layer
    [dy0,dy1,dy2] = system(t , y, PeT, Mx, Mp, Veff, gamma)
    """
    dy0=y[1]
    dy1=P/Mx*(V*Mp+(V+y[2])*y[1]/(2*y[0]*G)*(1-Mp*t-Mx*y[0])-Mx*V*y[1])
    dy2=-(V+y[2])*y[1]/(2*y[0])
    
    return [ dy0 , dy1 , dy2 ]



def systemTest(t, y):
    return -t**2+2*t

def systemTest2(t, y, a):
    return a

def systemTest3(t , y):
    return [y[1]+1, 2*y[0]]


