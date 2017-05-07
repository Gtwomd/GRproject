from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint

#FROM MATPLOTLIB WEBSITE
################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')
################################

def geo(y, s):
	theta = y[0]
	phi = y[1]
	dtheta = y[2]
	dphi = y[3]

	dyds = np.zeros_like(y)
	dyds[0] = dtheta
	dyds[1] = dphi
	dyds[2] = np.sin(theta)*np.cos(theta)*dphi**2
	dyds[3] = -2 *(np.tan(theta)**-1)*dphi*dtheta
	return dyds

x = np.linspace(0,4*np.pi,10000)

y = odeint(geo, [np.pi/3,0,.4,1], x)

theta = y[:,0]
phi = y[:,1]

ax.plot(
	np.cos(phi)*np.sin(theta),
	np.sin(phi)*np.sin(theta),
	np.cos(theta)
	)
plt.show()