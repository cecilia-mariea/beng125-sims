import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from scipy.integrate import solve_ivp
from itertools import product

BASE_PARAMETERS = {
    "kp" : 0.8, # p53 basal synthesis
    "kp_p53" : 6, # autoregulation of p53
    "dp" : 0.1, # natural degredation of p53
    "Kp" : 2, # half saturation constant for TP53
    "Km" : 1,  # half saturation constant for MDM2 mRNA
    "lam" : 3, # ubiquitination of p53
    "n" : 4, #  cooperativity of TP53 binding
    "m": 2, # cooperativity of MDM2 mRNA binding
    "km": 0.5, # MDM2 basal synthesis
    "km_p53" : 1.5, # activation of MDM2
    "dm" : 0.2 # natural degredation of MDM2
} 

def hill(p, n, K):
    p = np.asarray(p, dtype=float)
    return np.power(p, n) / (np.power(K, n) + np.power(p, n))

def dpdt(p, M, params):

    kp = params["kp"]
    kp_p53 = params["kp_p53"]
    dp = params["dp"]
    Kp = params["Kp"]
    lam = params["lam"]
    n = params["n"]

    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * M) * p)

def dMdt(p, M, params):

    Km = params["Km"]
    m = params["m"]
    km = params["km"]
    km_p53 = params["km_p53"]
    dm = params["dm"]

    return (km + km_p53 * hill(p, m, Km)) - (dm * M)

def p53_MDM2(t, vars, params):

    p, M = vars

    return (dpdt(p, M, params), dpdt(p, M, params))

def solve_sys(t_span, t_eval, init_cond, params):
    return solve_ivp(
        fun = lambda t, vars: p53_MDM2(t, vars, params), 
        t_span = t_span,
        y0 = init_cond, # should be array-like
        t_eval = t_eval,
        rtol=1e-8, # relative tolerance
        atol=1e-10 # global tolerance
    )

def t_sim(t0, tf, init_cond, params, init_range=False):

    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    
    if init_range:

        num_plots = len(init_cond)
        n = ceil(np.sqrt(num_plots))
        fig, axes = plt.subplots(n, n, figsize=(20,20))
        axes = axes.flatten()

        for ax in axes[num_plots:]: ax.axis('off')

        for curr_cond, ax in zip(init_cond, axes):

            sol = solve_sys(t_span, t_eval, curr_cond, params)

            ax.plot(sol.t, sol.y[0], label=f"p53 concentration, p53(t=0)={curr_cond[0]} uM")
            ax.plot(sol.t, sol.y[1], label=f"MDM2 concentration, MDM2(t=0)={curr_cond[1]} uM")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Species Concentration (uM)")
            ax.grid(True, alpha=0.3)
        
    else:

        plt.figure(figsize=(8,5))

        sol = solve_sys(t_span, t_eval, init_cond, params)

        plt.plot(sol.t, sol.y[0], label="p53 concentration")
        plt.plot(sol.t, sol.y[1], label="MDM2 concentration")

        plt.grid(True, alpha=0.3)
        plt.ylabel("Species Concentration (uM)")
        plt.xlabel("Time (min)")

    plt.legend()
    plt.tight_layout()
    plt.show()

def pairwise(p0, M0): 
    return tuple(product(p0, M0))

if __name__ == "__main__":

    # if many init cond, pass a itera of init cond pairs (p,M)

    p0 = np.linspace(1,5,5)
    M0 = np.linspace(1,5,5)
    init_cond = pairwise(p0,M0) 
    # init_cond = (1,2)
    t0 = 0
    tf = 15

    t_sim(t0, tf, init_cond, BASE_PARAMETERS, init_range=True)