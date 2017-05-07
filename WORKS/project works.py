from __future__ import division
import gravipy as gp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import lambdify, N, Function, Derivative, solve
from sympy.abc import a,b,c,d,e,f,g,h,i,j,k
from sympy.utilities.lambdify import lambdify, implemented_function

print '\n\n'

# m_earth = 6e24 #kilograms
# m = m_earth
# rs = 8.87e-3
rs = 2.95e3


t,r,theta,phi,tau = gp.symbols('t r \\theta \phi \\tau')
coordinates = gp.Coordinates('\chi', [t, r, theta, phi])
# Metric2D = gp.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, gp.sin(theta)*r**2)
Metric2D = gp.diag(-(1-rs/r), 1/(1-rs/r), r**2, gp.sin(theta)*r**2)

gtensor = gp.MetricTensor('g', coordinates, Metric2D)
w = gp.Geodesic('w', gtensor, tau)


d2t = lambda a,b,c,d,e,f,g,h: (
	solve(w(1),Derivative(t(tau),tau,tau))[0].subs(
		{
		t(tau):a,
		r(tau):b,
		theta(tau):c,
		phi(tau):d,
		Derivative(t(tau),tau):e,
		Derivative(r(tau),tau):f,
		Derivative(theta(tau),tau):g,
		Derivative(phi(tau),tau):h,
		}
		)
	)
d2r = lambda a,b,c,d,e,f,g,h: (
	solve(w(2),Derivative(r(tau),tau,tau))[0].subs(
		{
		t(tau):a,
		r(tau):b,
		theta(tau):c,
		phi(tau):d,
		Derivative(t(tau),tau):e,
		Derivative(r(tau),tau):f,
		Derivative(theta(tau),tau):g,
		Derivative(phi(tau),tau):h,
		}
		)
	)
d2theta = lambda a,b,c,d,e,f,g,h: (
	solve(w(3),Derivative(theta(tau),tau,tau))[0].subs(
		{
		t(tau):a,
		r(tau):b,
		theta(tau):c,
		phi(tau):d,
		Derivative(t(tau),tau):e,
		Derivative(r(tau),tau):f,
		Derivative(theta(tau),tau):g,
		Derivative(phi(tau),tau):h,
		}
		)
	)
d2phi = lambda a,b,c,d,e,f,g,h: (
	solve(w(4),Derivative(phi(tau),tau,tau))[0].subs(
		{
		t(tau):a,
		r(tau):b,
		theta(tau):c,
		phi(tau):d,
		Derivative(t(tau),tau):e,
		Derivative(r(tau),tau):f,
		Derivative(theta(tau),tau):g,
		Derivative(phi(tau),tau):h,
		}
		)
	)



# d2t = implemented_function(Function('d2t'), d2t)
d2t = lambdify((a,b,c,d,e,f,g,h), d2t(a,b,c,d,e,f,g,h), modules='numpy')

# d2r = implemented_function(Function('d2r'), d2r)
d2r = lambdify((a,b,c,d,e,f,g,h), d2r(a,b,c,d,e,f,g,h), modules='numpy')

# d2theta = implemented_function(Function('d2theta'), d2theta)
d2theta = lambdify((a,b,c,d,e,f,g,h), d2theta(a,b,c,d,e,f,g,h), modules='numpy')

# d2phi = implemented_function(Function('d2phi'), d2phi)
d2phi = lambdify((a,b,c,d,e,f,g,h), d2phi(a,b,c,d,e,f,g,h), modules='numpy')

# x = (1,2,3,4,5,6,7,8)

# print float(d2t(*x))
# print float(d2r(*x))
# print float(d2phi(*x))
# print float(d2theta(*x))


def geo(x, tau):
    '''
    x[0] = t
    x[1] = r
    x[2] = theta
    x[3] = phi
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



#####################################################
r_earth = 1.5e11 #meters
t_earth = 3e7 * 3e8 #seconds
drdtau_earth = 9e-2 #meters per second
dphidtau_earth = 2e-8 #radians per secong
#####################################################
t0 = 0
r0 = r_earth
theta0 = np.pi/2
phi0 = 0
# dt0 = 200
dt0 = 1
dr0 = 0
dtheta0 = 0
# dphi0 = 15e-17
dphi0 = 7e-16
#####################################################
x0_ = np.array([t0,r0,theta0,phi0,dt0,dr0,dtheta0,dphi0])
#####################################################
tau_ = np.linspace(0,2*t_earth, 10000)
x_ = odeint(geo, x0_, tau_)
#####################################################
t_ = x_[:,0]
r_ = x_[:,1]
theta_ = x_[:,2]
phi_ = x_[:,3]
#####################################################
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#MPL WEBSITE
#####################################################
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = r_earth* np.outer(np.cos(u), np.sin(v))
y = r_earth* np.outer(np.sin(u), np.sin(v))
z = r_earth* np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')
ax.set_xlabel('x')
ax.set_ylabel('y')
#####################################################
ax.plot(
	r_*np.cos(phi_)*np.sin(theta_),
	r_*np.sin(phi_)*np.sin(theta_),
	r_*np.cos(theta_)
	)

plt.show()

