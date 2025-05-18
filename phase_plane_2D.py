# imports
import numpy as np
from numpy.linalg import eigvals 
import matplotlib.pyplot as plt

from scipy.differentiate import jacobian
from scipy.integrate import odeint
from scipy.optimize import fsolve

from set_params_eqns import EXP_PARAMETERS_2D, dMdt, p53_MDM2 
from set_params_eqns import dpdt_2D as dpdt

def calc_vec_field():
    pass

def calc_nullclines():
    pass

# Find fixed points
def calc_fixed_pts(p0, m0):

    fixed_pts = []

    for p in p0: 
        for m in m0:
            initial_guess = np.array([p,m])
            sol = fsolve(p53_MDM2, initial_guess, args=(EXP_PARAMETERS_2D))

            if np.allclose(
                [dpdt(sol[0], sol[1], EXP_PARAMETERS_2D),
                 dMdt(sol[0], sol[1], EXP_PARAMETERS_2D)],
                [0,0]):
                fixed_pts.append(tuple(sol))

    return fixed_pts

def calc_jacobian():
    pass

def classify_fixed_pts(pts):
    pass

def plot_phase_plane():
    pass

def simulate_trajectories():
    pass

if __name__ == '__main__':
    p0 = (0,) 
    m0 = (0,)
    pts = calc_fixed_pts(p0, m0)