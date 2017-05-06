import numpy as np
from scipy.integrate import odeint

def time(x, tau):
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
	dxdtau[4] = (2*M*x[5]*x[4])/((2*M-x[1])*x[1])
	dxdtau[5] = 
	dxdtau[6] = 
	dxdtau[7] = 

