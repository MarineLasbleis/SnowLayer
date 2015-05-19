import numpy as np
from scipy.integrate import ode
import systemdiff


def resolutionSystem3(y0, t0, tmax, dt, PeT=1.e6, MX=0.17, MP=0.016, Veff=1.e-7, gamma=7.2,Atol=1.e-10,Rtol=1.e-12):
    """ Integration RK (de type Dormand Prince 4-5) of the equations system for an horizontal F-layer.
    resolutionSystem3(y0, t0, tmax, dt, *arg)
    Default parameters for the others arguments are :  PeT=1.e6, MX=0.17, MP=0.016, Veff=1.e7, gamma=7.2,Atol=1.e-10,Rtol=1.e-12
    """
    
    sol = ode(systemdiff.system).set_integrator('dopri5',atol=Atol,rtol=Rtol,nsteps=1e7)
    sol.set_initial_value(y0, t0)
    sol.set_f_params(PeT, MX, MP, Veff, gamma)

    
    y1 , y2 , y3 ,  t = y0[0]*np.ones(1), y0[1]*np.ones(1), y0[2]*np.ones(1), t0*np.ones(1)
    

    if dt>0 :
        while sol.successful() and sol.t < tmax:
            sol.integrate(sol.t+dt)
        # print("%g %g" % (r.t, r.y))
            y1 = np.append(y1, sol.y[0])
            y2 = np.append(y2, sol.y[1])
            y3 = np.append(y3, sol.y[2])
            t = np.append(t, sol.t)
    else:
        while sol.successful() and sol.t > tmax:
            sol.integrate(sol.t+dt)
        # print("%g %g" % (r.t, r.y))
            y1 = np.append(y1, sol.y[0])
            y2 = np.append(y2, sol.y[1])
            y3 = np.append(y3, sol.y[2])
            t = np.append(t, sol.t)
            
    return y1 , y2, y3, t
