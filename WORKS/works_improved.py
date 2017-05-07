from __future__ import division
import gravipy as gp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import lambdify, N, Function, Derivative, solve
from sympy.abc import a,b,c,d,e,f,g,h,i,j,k
from sympy.utilities.lambdify import lambdify, implemented_function

rs = 2.95e3 #schwarzchild radius of sun

t,r,theta,phi,tau = gp.symbols('t r \\theta \phi \\tau')
coordinates = gp.Coordinates('\chi', [t, r, theta, phi])
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

d2t = lambdify((a,b,c,d,e,f,g,h), d2t(a,b,c,d,e,f,g,h), modules='numpy')
d2r = lambdify((a,b,c,d,e,f,g,h), d2r(a,b,c,d,e,f,g,h), modules='numpy')
d2theta = lambdify((a,b,c,d,e,f,g,h), d2theta(a,b,c,d,e,f,g,h), modules='numpy')
d2phi = lambdify((a,b,c,d,e,f,g,h), d2phi(a,b,c,d,e,f,g,h), modules='numpy')



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
c = 3e8
r_earth = 1.5e11 #meters
t_earth = 3e7 #seconds
drdtau_earth = 9e-2 #meters per second
dphidtau_earth = 2e-8 #radians per secong
#####################################################
t0 = 0
r0 = r_earth
theta0 = np.pi/2
phi0 = 0
dt0 = 1
dr0 = 0
dtheta0 = 0
dphi0 = 2e-7
#####################################################
x0_ = np.array([c*t0,r0,theta0,phi0,c*dt0/c,dr0/c,dtheta0/c,dphi0/c])
#####################################################
tau_ = np.linspace(0,2*c*t_earth, 10000)
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

#FROM MATPLOTLIB WEBSITE PLOT A SPHERE
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

