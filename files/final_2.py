import gravipy
import numpy as np
import matplotlib.pyplot
import sympy
from scipy.integrate import odeint
from sympy import lambdify, N, Function
from sympy.abc import a,b,c,d,e,f,g,h,i,j,k
from sympy.utilities.lambdify import lambdify, implemented_function
import sys

print '\n\n\n\n\n\n'

t,r,theta,phi,tau = gravipy.symbols('t r \\theta \phi \\tau')
c_ = gravipy.Coordinates('\chi', [t, r, theta, phi])

m = 5

Metric2D = gravipy.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, r**2)

g_ = gravipy.MetricTensor('g', c_, Metric2D)

w = gravipy.Geodesic('w', g_, tau)

d2t = lambda a,b,c,d,e,f,g,h: sympy.solve(w(1),sympy.Derivative(t(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(theta(tau),tau):c,sympy.Derivative(phi(tau),tau):d,t(tau):e,r(tau):f,theta(tau):g,phi(tau):h})
d2r = lambda a,b,c,d,e,f,g,h: sympy.solve(w(2),sympy.Derivative(r(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(theta(tau),tau):c,sympy.Derivative(phi(tau),tau):d,t(tau):e,r(tau):f,theta(tau):g,phi(tau):h})
d2theta = lambda a,b,c,d,e,f,g,h: sympy.solve(w(3),sympy.Derivative(theta(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(theta(tau),tau):c,sympy.Derivative(phi(tau),tau):d,t(tau):e,r(tau):f,theta(tau):g,phi(tau):h})
d2phi = lambda a,b,c,d,e,f,g,h: sympy.solve(w(4),sympy.Derivative(phi(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(theta(tau),tau):c,sympy.Derivative(phi(tau),tau):d,t(tau):e,r(tau):f,theta(tau):g,phi(tau):h})


d2t = implemented_function(Function('d2t'), d2t)
d2t = lambdify((a,b,c,d,e,f,g,h), d2t(a,b,c,d,e,f,g,h))

d2r = implemented_function(Function('d2r'), d2r)
d2r = lambdify((a,b,c,d,e,f,g,h), d2r(a,b,c,d,e,f,g,h))

d2theta = implemented_function(Function('d2theta'), d2theta)
d2theta = lambdify((a,b,c,d,e,f,g,h), d2theta(a,b,c,d,e,f,g,h))

d2phi = implemented_function(Function('d2phi'), d2phi)
d2phi = lambdify((a,b,c,d,e,f,g,h), d2phi(a,b,c,d,e,f,g,h))
print float(d2t(1,2,3,4,5,6,7,8))
print float(d2r(1,2,3,4,5,6,7,8))
print float(d2phi(1,2,3,4,5,6,7,8))
print float(d2theta(1,2,3,4,5,6,7,8))

def geo(x, tau):
    '''
    x[0] = t
    x[1] = r
    x[2] = phi
    x[3] = theta
    x[4] = dt dtau
    x[5] = dr dtau
    x[6] = dtheta dtau
    x[7] = dphi dtau
    '''
    dxdtau = np.zeros_like(x)

    dxdtau[0] = x[4]
    dxdtau[1] = x[5]
    dxdtau[2] = x[6]
    dxdtau[3] = x[7]
    dxdtau[4] = float(d2t(*x))
    dxdtau[5] = float(d2r(*x))
    dxdtau[6] = float(d2theta(*x))
    dxdtau[7] = float(d2phi(*x))
    return dxdtau

tau = np.linspace(0,1,2)
x0_ = np.array([1,2,0,4,5,6,0,6])
x_ = odeint(geo, x0_, tau)

print x_

