# imports
import numpy as np
from numpy.linalg import eigvals 
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
from scipy.differentiate import jacobian

# custom imports
from set_params_eqns import  dMdt, p53_MDM2 
from set_params_eqns import dpdt_2D as dpdt, PARAMETERS_2D as params

# find fixed points
def calc_fixed_pts(p0, m0):

    fixed_pts = []

    for p in p0: 
        for m in m0:

            initial_guess = np.array([p,m])
            sol = fsolve(p53_MDM2, initial_guess, args=(params))

            p_star = sol[0]
            M_star = sol[1]

            # check if solver was successful
            if np.allclose(
                [dpdt(p_star, M_star, params),
                 dMdt(p_star, M_star, params)],
                [0,0]):

                # round solutions to avoid non-unique solutions from overflow
                sol_rounded = (round(p_star,4), round(M_star,4))
                fixed_pts.append(sol_rounded)

    return set(fixed_pts)

# TODO
def classify_fixed_pts(pts):

    def calc_jacobian():
        pass

    pass

# TODO
def plot_nullclines():
    pass

if __name__=="__main__":
    print("plot nullclines script called directly")