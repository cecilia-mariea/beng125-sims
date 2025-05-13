# imports
import numpy as np
import matplotlib.pyplot as plt

from scipy.differentiate import jacobian
from scipy.integrate import odeint
from scipy.optimize import fsolve

from diffeq_2D import BASE_PARAMETERS, p53_MDM2, solve_sys,dpdt, dMdt

# Find fixed points
def find_fixed_pts(p0, m0):

    fixed_pts = []

    for p in p0: 
        for m in m0:
            initial_guess = np.array([p,m])
            sol = fsolve(p53_MDM2, initial_guess)

            if np.allclose(
                [dpdt(sol[0], sol[1], BASE_PARAMETERS),
                 dMdt(sol[0], sol[1], BASE_PARAMETERS)],
                [0,0]):
                fixed_pts.append(tuple(sol))

    return fixed_pts

# Plot Phase Plane
def plot_phase_plane(nullcline=False):
    if nullcline:
        pass

if __name__ == '__main__':
    pass
