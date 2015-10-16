# Time-stamp: <2015-10-13 17:11:08 marine>
## Project : Snow in the F - layer
## Subproject : computation of the equilibrium state and stability of it
## Author : Marine Lasbleis


## Tests


########### Import external libraries + own files : OK
########### Use of parameters in function f : OK !

# external libraries
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#files
import systemdiff

y0, t0 = 0.0 , 0
y , t = y0*np.ones(1) , t0*np.ones(1)
b = 1

sol = ode(systemdiff.systemTest2).set_integrator('dopri5')
sol.set_initial_value(y0, t0)
sol.set_f_params(b)
t1 = 10
dt = 1
i = 1

while sol.successful() and sol.t < t1:
    sol.integrate(sol.t+dt)
   # print("%g %g" % (r.t, r.y))
    y = np.append(y, sol.y)
    t = np.append(t, sol.t)

print "=============="
print "=============="
print "Test 1 : using system of equations systemTest2, with y[t=0]=0."
print "With the given parameter b=",b","
print "the solution should be y=bx"
print "=============="
print "=============="
    
plt.plot(t, y)
#plt.plot(t,np.exp(t))
plt.xlabel('t')
plt.ylabel('y')
plt.show()


##########
##########

#Test2 : y is an array of length 2  : OK !

# external libraries
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#files
import systemdiff

y0, t0 = [1.0 , 2.0], 0
y1 , y2 , t = y0[0]*np.ones(1), y0[1]*np.ones(1), t0*np.ones(1)

sol = ode(systemdiff.systemTest3).set_integrator('dopri5')
sol.set_initial_value(y0, t0)

t1 = 1
dt = 0.01
i = 1

while sol.successful() and sol.t < t1:
    sol.integrate(sol.t+dt)
   # print("%g %g" % (r.t, r.y))
    y1 = np.append(y1, sol.y[0])
    y2 = np.append(y2, sol.y[1])
    t = np.append(t, sol.t)

print "=============="
print "=============="
print "Test2 : same, with array of length 2!"
print "it solves the system of equation systemTest3, with y[0]=[1,2]"
print "solution is y1(t)=exp(2t), y2(t)=2exp(2t)"
print ""+" is exp(2t)"
print "=============="
print "=============="

    
plt.plot(t, y1, t, y2)
plt.plot(t,np.exp(2.*t),'+')
plt.xlabel('t')
plt.ylabel('y')
plt.show()





##########
##########

#Test3 : y is an array of length 3

# external libraries
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt

#files
import systemdiff

y0, t0 = [3.0 , 1.0, 1.0], 0
y1 , y2 , y3 ,  t = y0[0]*np.ones(1), y0[1]*np.ones(1), y0[2]*np.ones(1), t0*np.ones(1)

sol = ode(systemdiff.systemTest4).set_integrator('dopri5')
sol.set_initial_value(y0, t0)

t1 =1
dt = 0.01
i = 1

while sol.successful() and sol.t < t1:
    sol.integrate(sol.t+dt)
   # print("%g %g" % (r.t, r.y))
    y1 = np.append(y1, sol.y[0])
    y2 = np.append(y2, sol.y[1])
    y3 = np.append(y3, sol.y[2])
    t = np.append(t, sol.t)


print "=============="
print "=============="
print "Test3 : same, with array of length 3!"
print "=============="
print "=============="
    
plt.plot(t, y1, t, y2, t, y3)
#plt.plot(t,np.exp(t))
plt.xlabel('t')
plt.ylabel('y')
plt.show()
