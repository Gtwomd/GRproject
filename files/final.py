import gravipy
import numpy as np
import matplotlib.pyplot
import sympy
from scipy.integrate import odeint
from sympy import lambdify
from sympy.abc import a,b,c,d
import sys
t,r,theta,phi,tau = gravipy.symbols('t r \\theta \phi \\tau')
c = gravipy.Coordinates('\chi', [t, r, theta, phi])

m = 5

Metric2D = gravipy.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, r**2)

g = gravipy.MetricTensor('g', c, Metric2D)

w = gravipy.Geodesic('w', g, tau)


d2t = sympy.solve(w(1),sympy.Derivative(t(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(phi(tau),tau):c,sympy.Derivative(theta(tau),tau):d})
d2r = sympy.solve(w(2),sympy.Derivative(r(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(phi(tau),tau):c,sympy.Derivative(theta(tau),tau):d})
d2theta = sympy.solve(w(3),sympy.Derivative(theta(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(phi(tau),tau):c,sympy.Derivative(theta(tau),tau):d})
d2phi = sympy.solve(w(4),sympy.Derivative(phi(tau),tau,tau))[0].subs({sympy.Derivative(t(tau),tau):a,sympy.Derivative(r(tau),tau):b,sympy.Derivative(phi(tau),tau):c,sympy.Derivative(theta(tau),tau):d})

d2tf = lambdify((a,b,c,d,t(tau),r(tau),theta(tau),phi(tau)),d2t,modules="numpy")
d2rf = lambdify((a,b,c,d,t(tau),r(tau),theta(tau),phi(tau)),d2r,modules="numpy")
d2thetaf = lambdify((a,b,c,d,t(tau),r(tau),theta(tau),phi(tau)),d2theta,modules="numpy")
d2phif = lambdify((a,b,c,d,t(tau),r(tau),theta(tau),phi(tau)),d2phi,modules="numpy")


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
    dxdtau[4] = d2tf(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7])
    dxdtau[5] = d2rf(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7])
    dxdtau[6] = d2thetaf(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7])
    dxdtau[7] = d2phif(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7])
    return dxdtau

tau = np.arange(0,10,.1)
x0_ = np.array([1,2,0,4,5,6,0,6])
x_ = odeint(geo, x0_, tau)

print x_

