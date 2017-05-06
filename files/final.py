import gravipy
import numpy as np
import matplotlib.pyplot
import sympy
from scipy.integrate import odeint
from sympy import lambdify
from sympy.abc import a,b
import sys
t,r,theta,phi,tau = gravipy.symbols('t r \\theta \phi \\tau')
c = gravipy.Coordinates('\chi', [t, r, theta, phi])

m = 5

Metric2D = gravipy.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, r**2)

g = gravipy.MetricTensor('g', c, Metric2D)

w = gravipy.Geodesic('w', g, tau)

d2t = lambda x: sympy.solve(w(1),sympy.Derivative(t(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):x[0],sympy.Derivative(r(tau),tau):x[1],sympy.Derivative(theta(tau),tau):x[2],sympy.Derivative(phi(tau),tau):x[3],t(tau):x[4],r(tau):x[5],theta(tau):x[6],phi(tau):x[7]})
d2r = lambda x: sympy.solve(w(2),sympy.Derivative(r(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):x[0],sympy.Derivative(r(tau),tau):x[1],sympy.Derivative(theta(tau),tau):x[2],sympy.Derivative(phi(tau),tau):x[3],t(tau):x[4],r(tau):x[5],theta(tau):x[6],phi(tau):x[7]})
d2theta = lambda x: sympy.solve(w(3),sympy.Derivative(theta(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):x[0],sympy.Derivative(r(tau),tau):x[1],sympy.Derivative(theta(tau),tau):x[2],sympy.Derivative(phi(tau),tau):x[3],t(tau):x[4],r(tau):x[5],theta(tau):x[6],phi(tau):x[7]})
d2phi = lambda x: sympy.solve(w(4),sympy.Derivative(phi(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):x[0],sympy.Derivative(r(tau),tau):x[1],sympy.Derivative(theta(tau),tau):x[2],sympy.Derivative(phi(tau),tau):x[3],t(tau):x[4],r(tau):x[5],theta(tau):x[6],phi(tau):x[7]})

print d2t([1,2,3,4,5,6,7,8])
print d2r([1,2,3,4,5,6,7,8])
print d2phi([1,2,3,4,5,6,7,8])
print d2theta([1,2,3,4,5,6,7,8])

def geo(x, tau):
    '''
    x[0] = t
    x[1] = r
    x[2] = phi
    x[3] = theta
    x[4] = dt dtau
    x[5] = dr dtau
    x[6] = dphi dtau
    x[7] = dtheta dtau
    '''
    dxdtau = np.zeros_like(x)
    dxdtau[0] = x[4]
    dxdtau[1] = x[5]
    dxdtau[2] = x[6]
    dxdtau[3] = x[7]
    dxdtau[4] = d2t(x)
    dxdtau[5] = d2r(x)
    dxdtau[6] = d2theta(x)
    dxdtau[7] = d2phi(x)
    return dxdtau

tau = [.1]#np.linspace(0,1,2)
x0_ = np.array([1,2,0,4,5,6,0,6])
x_ = odeint(geo, x0_, float(sys.argv[1]))

print x_

