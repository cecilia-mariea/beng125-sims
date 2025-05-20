# imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# custom imports
from set_params_eqns import  dMdt, p53_MDM2, calc_nullcline
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
                sol_rounded = round(p_star,4), round(M_star,4)
                fixed_pts.append(sol_rounded)

    return set(fixed_pts)

def jacobian(pt, h=1e-6): 

    p, M = pt

    partial_M_partial_M = (dMdt((p, M+h), **params) - dMdt(p, M-h **params)) / 2*h
    partial_M_partial_p = (dMdt((p+h, M), **params) - dMdt((p-h, M), **params)) / 2*h
    partial_p_partial_M = (dpdt((p,M+h),**params) - dpdt((p,M-h), **params)) / 2*h
    partial_p_partial_p =  (dpdt((p+h, M), **params) - dpdt((p-h, M), **params)) / 2*h

    return np.array([[partial_p_partial_p, partial_p_partial_M],
            [partial_M_partial_p, partial_M_partial_M]])

def classify_fixed_pts(pts):

    classified_pts = np.zeros_like(pts)

    def classify(eigenvals):

        if np.all(np.real(eigenvals)) > 0:
            return 'unstable'
        elif np.all(np.real(eigenvals)) < 0:
            return 'stable'
        elif (np.any(np.real(eigenvals)) < 0 ) and (np.any(np.real(eigenvals > 0))):
            return 'saddle node'
        elif (np.any(np.img(eigenvals)) != 0):
            return 'oscillatory'
        else:
            return 'linear analysis fails'

    for i, pt in enumerate(pts):

        J = jacobian(pt)
        eigenvals = np.linalg.eigvals(J)

        classified_pts[i] = pt, classify(eigenvals)

    return classified_pts

def plot_nullclines(p):

    P, p_nullcline_M, M_nullcline_M = calc_nullcline(p, **params)

    plt.figure(figsize=(10,10))

    plt.plot(P, p_nullcline_M, label="P53 nullcline", color="r", linewidth=3)
    plt.plot(P, M_nullcline_M, label="MDM2 nullcline", color="b", linewidth=3)

    plt.xlabel("dpDt")
    plt.ylabel("dMdt")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.show()

if __name__=="__main__":
    print("plot nullclines script called directly")