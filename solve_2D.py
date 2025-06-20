import numpy as np
import matplotlib.pyplot as plt
from math import ceil
from scipy.integrate import solve_ivp
from itertools import product

from set_params_eqns import p53_MDM2_dt
from set_params_eqns import PARAMETERS_2D as params

def solve_sys(t_span, t_eval, init_cond):
    return solve_ivp(
        fun = lambda t, vars: p53_MDM2_dt(t, vars, params), 
        t_span = t_span,
        y0 = init_cond, # should be array-like
        t_eval = t_eval,
        dense_output=True,
        rtol = 1e-8, # relative tolerance
        atol = 1e-10 # global tolerance
    )

def t_sim(t0, tf, init_cond, init_range=False):

    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)

    sols = []
    
    if init_range:

        num_plots = len(init_cond)
        n = ceil(np.sqrt(num_plots))
        fig, axes = plt.subplots(n, n, figsize=(20,20))
        axes = axes.flatten()

        for ax in axes[num_plots:]: ax.axis('off')

        for curr_cond, ax in zip(init_cond, axes):

            sol = solve_sys(t_span, t_eval, curr_cond)
            sols.append(sol)

            ax.plot(sol.t, sol.y[0], label=f"p53 concentration, p53(t=0)={curr_cond[0]} uM")
            ax.plot(sol.t, sol.y[1], label=f"MDM2 concentration, MDM2(t=0)={curr_cond[1]} uM")
            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Species Concentration (uM)")
            ax.grid(True, alpha=0.3)
        
    else:

        plt.figure(figsize=(8,5))

        sol = solve_sys(t_span, t_eval, init_cond)
        sols.append(sol)

        plt.plot(sol.t, sol.y[0], label="p53 concentration")
        plt.plot(sol.t, sol.y[1], label="MDM2 concentration")

        plt.grid(True, alpha=0.3)
        plt.ylabel("Species Concentration (uM)")
        plt.xlabel("Time (min)")
        plt.title(f"Time series of Species Concentration with Initial Conditions {init_cond}")


    plt.legend()
    plt.tight_layout()
    plt.show()
    
    return sol

pairwise = lambda p0, M0: tuple(product(p0, M0))

if __name__ == "__main__":

    print("2D solve script called explicitly")
