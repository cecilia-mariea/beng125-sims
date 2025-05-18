# imports
import numpy as np
from numpy.linalg import eigvals 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.differentiate import jacobian

# module imports
from set_params_eqns import  dMdt, p53_MDM2 
from set_params_eqns import dpdt_2D as dpdt, EXP_PARAMETERS_2D as params

def calc_vec_field(p_range, M_range, grid_size):

    p = np.linspace(p_range[0], p_range[1], grid_size)
    m = np.linspace(M_range[0], M_range[1], grid_size)
    P, M = np.meshgrid(p,m)

    dpdt_spec_sol = np.zeros_like(P)
    dMdt_spec_sol = np.zeros_like(M)

    for i in range(grid_size):
        for j in range(grid_size):

            vars = (P[i,j], M[i,j])

            dpdt_spec_sol[i,j] = dpdt(vars, **params)
            dMdt_spec_sol[i,j] = dMdt(vars, **params)
    
    magnitudes = np.sqrt(dpdt_spec_sol**2 + dMdt_spec_sol**2)

    return P, M, dpdt_spec_sol, dMdt_spec_sol, magnitudes

def calc_nullclines():
    pass

# Find fixed points
def calc_fixed_pts(p0, m0):

    fixed_pts = []

    for p in p0: 
        for m in m0:
            initial_guess = np.array([p,m])
            sol = fsolve(p53_MDM2, initial_guess, args=(params))

            p_star = sol[0]
            M_star = sol[1]

            if np.allclose(
                [dpdt(p_star, M_star, params),
                 dMdt(p_star, M_star, params)],
                [0,0]):

                sol_rounded = (round(p_star,4), round(M_star,4))
                fixed_pts.append(sol_rounded)

    return set(fixed_pts)

def calc_jacobian():
    pass

def classify_fixed_pts(pts):
    pass

def plot_phase_plane():
    pass

def simulate_trajectories():
    pass

if __name__ == '__main__':
    print("2D phase plane script was explicitly called")