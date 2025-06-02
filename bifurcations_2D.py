# imports
import numpy as np
from numpy.linalg import eigvals
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.differentiate import jacobian

# custom imports
from set_params_eqns import dMdt, p53_MDM2
from set_params_eqns import dpdt_2D as dpdt, PARAMETERS_2D as PARAMS

GUESSES = np.array([[x, y] for x in [0, 5, 10, 50] for y in [0, 5, 10, 50]])

def calc_fixed_pts(a, b, a_val, b_val, guesses=GUESSES):

    params = PARAMS.copy()
    params[a] = a_val
    params[b] = b_val

    pts = []

    for guess in guesses:
        sol = fsolve(p53_MDM2, guess, args=(params))

        p_star = sol[0]
        M_star = sol[1]
        vec = (p_star, M_star)

        # check if solver was successful
        if np.allclose(
            [dpdt(vec, **params),
                dMdt(vec, **params)],
            [0,0]):

            if not any(np.allclose(pt, vec, rtol=1e-5, atol=1e-8) for pt in pts):

                pts.append(vec)

    return pts

def classify_fixed_pt(pt, a, b, a_val, b_val, simple=False, spirals=True):

    def classify(eig):
        if np.all(np.iscomplex(eig)):
            if np.all(eig.real < 0):
                return "stable spiral"
            elif np.all(eig.real > 0):
                return "unstable spiral"
            elif np.allclose(eig.real, 0):
                return "center"
        else:
            if np.any(eig.real < 0) and np.any(eig.real > 0):
                return "saddle"
            elif np.all(eig.real < 0):
                return "stable node"
            elif np.all(eig.real > 0):
                return "unstable node"
        return "other"

    def simple_classify(eig):
        if np.all(np.iscomplex(eig)):
            return "complex"
        else:
            return "real"

    def spiral_classify(eig):
        if np.all(np.iscomplex(eig)):
            if np.all(eig.real < 0):
                return "stable"
            elif np.all(eig.real > 0):
                return "unstable"
            elif np.allclose(eig.real, 0):
                return "cycle"
    
    params = PARAMS.copy()
    params[a] = a_val
    params[b] = b_val

    eqn = lambda vars: (dpdt(vars, **params), dMdt(vars,**params))

    J = jacobian(
        eqn,
        pt,
        tolerances = {"rtol": 1e-8, "atol": 1e-10}
    )

    eigenvals = np.linalg.eigvals(J.df)

    if simple:
        return simple_classify(eigenvals)
    elif spirals:
        return spiral_classify(eigenvals)
    else:
        return classify(eigenvals)