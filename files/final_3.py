from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

r_earth = 150e9 #meters
t_earth = 3e7 * 3e8 #seconds
drdtau_earth = 9e-2 #meters per second
# dphidtau_earth = 7.3e-5 #radians per secong
dphidtau_earth = 1


rs = 8.87e-3
c=3e8

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
    dxdtau[4] = -(rs/(x[1]*(x[1]-rs)))*x[4]*x[5]
    # dxdtau[5] = -(
    # 	+(1/2)*(-rs/(x[1]*(x[1]-rs)))*(x[5])**2
    # 	-(x[1]-rs)*(x[6])**2
    # 	-(x[1]-rs)*(np.sin(x[2])**2)*(x[7])**2
    # 	+(1/2)*(c**2)*((rs*(x[1]-rs))/(x[1]**3))*(x[4])**2
    # 	)
    dxdtau[5] = -(
    	+(1/2)*(-rs/(x[1]*(x[1]-rs)))*(x[5])**2
    	-(x[1]-rs)*(x[7])**2
    	+(1/2)*(c**2)*((rs*(x[1]-rs))/(x[1]**3))*(x[4])**2
    	)
    # dxdtau[6] = -(
    # 	+(2/x[1])*x[6]*x[5]
    # 	-np.sin(x[2])*np.cos(x[2])*x[7]**2
    # 	)
    dxdtau[6] = 0
    # dxdtau[7] = -(
    # 	+(2/x[1])*x[7]*x[5] + 2*(np.tan(x[2])**-1)*x[6]*x[7]
    # 	)
    dxdtau[7] = -(2/x[1])*x[7]*x[5]
    return dxdtau

x0_ = np.array([
	0,
	r_earth,
	np.pi/2,
	0,
	1,
	0,
	0,
	dphidtau_earth
	])


tau_ = np.linspace(0,t_earth)
x_ = odeint(geo, x0_, tau_)

print x_

ax = plt.subplot(111, projection='polar')
ax.plot(x_[:,3], x_[:,1])
plt.show()

