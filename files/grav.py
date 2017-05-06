import gravipy
import numpy as np
import matplotlib.pyplot
import sympy
# gravipy.init_printing()

t, m,r,theta,phi,tau = gravipy.symbols('t M r \\theta \phi \\tau')

x = gravipy.Coordinates('\chi', [t, r, theta, phi])
Metric = gravipy.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, r**2*gravipy.sin(theta)**2)
Metric2D = gravipy.diag(-(1-2*m/r), 1/(1-2*m/r), r**2, r**2)
Metric_weak = gravipy.diag(-(1+2*m/r), (1+2*m/r), (1+2*m/r), (1+2*m/r))

g = gravipy.MetricTensor('g', x, Metric2D)
w = gravipy.Geodesic('w', g, tau)
print w(gravipy.All)