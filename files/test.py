from __future__ import division
import numpy as np
from sympy import tan, cos, sin, Matrix
from sympy.abc import u,v 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint



C =  Matrix([[(0, -tan(v)), (0,0)],[(sin(v)*cos(v),0),(0, 0)]])

def f(y,s,C,u,v):
    y0 = y[0] # u
    y1 = y[1] # u'
    y2 = y[2] # v
    y3 = y[3] # v'
    dy = np.zeros_like(y)
    dy[0] = y1
    dy[2] = y3

    C = C.subs({u:y0,v:y2}) # Evaluate C for u,v = (u0,v0)

    # dy[1] = -C[0,0][0]*dy[0]**2 -\
    #        2*C[0,0][1]*dy[0]*dy[2] -\
    #          C[0,1][1]*dy[2]**2
    dy[1] = 0
    dy[3] = -C[1,0][0]*dy[0]**2 -\
           2*C[1,0][1]*dy[0]*dy[2] -\
             C[1,1][1]*dy[2]**2
    return dy

def solve(C,u0,s0,s1,ds):
    s = np.arange(s0,s1+ds,ds)
    # The Christoffel symbol of 2nd kind, C, is a function of (u,v)
    
    return odeint(f,u0,s,args=(C,u,v))  # integration method : LSODA

y = solve(C,[0,1,np.pi/2,0],0,10,.01)

s = np.arange(0,10.01,.01)

print y[:,2]
print s

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(
    np.cos(y[:,2])*np.sin(y[:,0]),
    np.sin(y[:,2])*np.sin(y[:,0]),
    np.cos(y[:,0])
    )

# plt.plot(y[:,2])

plt.show()