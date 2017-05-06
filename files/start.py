import numpy as np
import matplotlib.pyplot as plt
import gravipy as gp
from scipy.integrate import odeint

def osc(y, t):
	'''
	y[0] = postion
	y[1] = velocity
	'''
	dydt = np.zeros_like(y)
	dydt[0] = y[1]
	dydt[1] = -10*y[0] 
	return dydt

x = np.linspace(0,100)
y = odeint(osc, [0,10], x)

plt.plot(x, y[:,0])
plt.show()

def coupled(y, t):
	dydt = np.zeros_like(y)
	dydt[0] = 5*y[0] + 3*y[1]
	dydt[1] = y[0] + 7*y[1]
	return dydt

x = np.linspace(0,1)
y = odeint(coupled, [5,1], x)

plt.plot(y[:,0], y[:,1])
plt.show()