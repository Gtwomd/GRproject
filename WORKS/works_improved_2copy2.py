from __future__ import division
import gravipy as gp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as colors
from scipy.integrate import odeint
from sympy import lambdify, N, Function, Derivative, solve
from sympy.abc import a,b,c,d,e,f,g,h,i,j,k
from sympy.utilities.lambdify import lambdify, implemented_function
import sys

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
#redefine c so as not to conflict with our import
c = 3e8
#####################################################
def trajectory(r_peri, v_peri, tau_):
	t0 = 0 # start at t=0
	r0 = r_peri
	theta0 = np.pi/2
	phi0 = 0
	dt0 = 1
	dr0 = 0 #no radial velocity at perihelion
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
	g_ = x_[:,4]

	return (t_, r_, theta_, phi_, g_)

#####################################################
#calculate over one earth year
days = 365
percent_year = days/365
earth_period = 3e7 #seconds
percent_year = (percent_year if len(sys.argv) <= 2 else float(sys.argv[2]))
tau_ = np.linspace(0,percent_year*c*earth_period, 10000*percent_year)
#####################################################
#format planet = (r_peri [in meters], v_peri [in meters per second])
mercury = (4.6e10, 5.898e4, 'Mercury','r') 
venus = (1.0748e11, 3.526e4, 'Venus','g')
earth = (1.471e11, 2.98e4, 'Earth','b')
mars = (2.067e11, 2.65e4, 'Mars','m')
planets_ = [mercury, venus, earth, mars]
#####################################################
#calculate trajectories
trajects = []
for planet in planets_:
	(t_, r_, theta_, phi_, g_) = trajectory(planet[0],planet[1], tau_)
	trajects.append([planet[2], t_, r_, theta_, phi_, g_, planet[3]])
#####################################################
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-2.5e11,2.5e11)
ax.set_ylim(-2.5e11,2.5e11)
ax.set_zlim(-2.5e11,2.5e11)

#FROM MATPLOTLIB WEBSITE PLOT A SPHERE
#plot sun
#####################################################
# Make data
sun_radius = 1.4e9
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = sun_radius * np.outer(np.cos(u), np.sin(v))
y = sun_radius * np.outer(np.sin(u), np.sin(v))
z = sun_radius * np.outer(np.ones(np.size(u)), np.cos(v))
# Plot the surface
ax.plot_surface(x, y, z, color='y')
ax.set_xlabel('meters')
ax.set_ylabel('meters')
ax.set_zlabel('meters')
#####################################################
trajec_data = []
for traject in trajects:
	(name, t_, r_, theta_, phi_, g_, color) = traject
	x = r_*np.cos(phi_)*np.sin(theta_)
	y = r_*np.sin(phi_)*np.sin(theta_)
	z = r_*np.cos(theta_)
	trajec_data.append([x,y,z,g_])
trajec_data = np.array(trajec_data)

#####################################################
def update_data(iteration, data, planets, labels, speed):
	iteration = np.minimum(iteration*speed,len(tau_)-1)
	g = "Gammas: "
	for num,planet in enumerate(planets):
		#set planet position
		planet.set_data(data[num,0,iteration], data[num,1,iteration])
		planet.set_3d_properties(data[num,2,iteration])
		# labels[num].set_data(data[num,0,iteration], data[num,1,iteration])
		# labels[num].set_3d_properties(data[num,2,iteration])
		#print gamma
		g+= str(planets_[num][2])+str(" : %.2e "%(data[num,3,iteration]-1.))
		#
		
	print g
	return planets
#####################################################
planet_draw = []
labels = []
for i in xrange(trajec_data.shape[0]):
	label = ax.text([],[],[], name)
	labels.append(label)
	j, = ax.plot([],[],[],marker="o")
	planet_draw.append(j)	

anim = animation.FuncAnimation(fig, update_data,100 , fargs=(trajec_data, planet_draw, labels, int(sys.argv[1])), interval=50, blit=True)

plt.show()

