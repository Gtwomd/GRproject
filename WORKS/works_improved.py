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

c = 3e8
#####################################################
def trajectory(r_peri, v_peri, tau_):
	t0 = 0
	r0 = r_peri
	theta0 = np.pi/2
	phi0 = 0
	dt0 = 1
	dr0 = 0
	dtheta0 = 0
	dphi0 = v_peri/r_peri
	

	x0_ = np.array([
		c*t0,
		r0,
		theta0,
		phi0,
		c*dt0/c,
		dr0/c,
		dtheta0/c,
		dphi0/c
		])

	x_ = odeint(geo, x0_, tau_)

	t_ = x_[:,0]
	r_ = x_[:,1]
	theta_ = x_[:,2]
	phi_ = x_[:,3]

	return (t_, r_, theta_, phi_)

#####################################################
earth_r_peri = 1.471e11 #meters
earth_v_peri = 2.98e4 #meters per second
earth_period = 3e7 #seconds
earth_tau = np.linspace(0,1*c*earth_period, 10000)

(t_, r_, theta_, phi_) = trajectory(earth_r_peri, earth_v_peri, earth_tau)

#####################################################
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(16,9))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-2e11,2e11)
ax.set_ylim(-2e11,2e11)
ax.set_zlim(-2e11,2e11)

#FROM MATPLOTLIB WEBSITE PLOT A SPHERE
#####################################################
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1e10 * np.outer(np.cos(u), np.sin(v))
y = 1e10 * np.outer(np.sin(u), np.sin(v))
z = 1e10 * np.outer(np.ones(np.size(u)), np.cos(v))
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

