# imports
import numpy as np
from scipy.optimize import fsolve
from scipy.differentiate import jacobian

# custom imports
from set_params_eqns import  dMdt, p53_MDM2
from set_params_eqns import dpdt_2D as dpdt, PARAMETERS_2D as params

GUESSES = np.array([[x, y] for x in [0, 5, 10, 50] for y in [0, 5, 10, 50]])

# find fixed points
def calc_fixed_pts(guesses=GUESSES):

    fixed_pts = []

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

            if not any(np.allclose(pt, vec, rtol=1e-5, atol=1e-8) for pt in fixed_pts):

                fixed_pts.append(vec)

    return fixed_pts

def classify_fixed_pt(pt):

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

    eqn = lambda vars: (dpdt(vars, **params), dMdt(vars, **params))

    J = jacobian(
        eqn,
        pt,
        tolerances = {"rtol": 1e-8, "atol": 1e-10}
    )

    eigenvals = np.linalg.eigvals(J.df)

    return classify(eigenvals)

if __name__=="__main__":
    print("fixed points script called directly")