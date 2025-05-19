import numpy as np
import matplotlib.pyplot as plt

from set_params_eqns import PARAMETERS_1D, dp

def dpdt_vs_p(p_vals, param_label="", param_vals=[]):

    plt.figure(figsize=(8,5))

    if bool(param_label):

        for val in param_vals: 

            params = PARAMETERS_1D.copy()
            params[param_label] = val

            dpdt_sol = np.array([dp(p, **params) for p in p_vals])

            plt.plot(p_vals, dpdt_sol, label=f"{param_label}={val}", alpha=0.7)
            plt.title(f'dp/dt versus p, Varying {param_label}')

    else: 

        dpdt_sol = np.array([dp(p, **PARAMETERS_1D) for p in p_vals] )

        plt.plot(p_vals, dpdt_sol)
        plt.title('dp/dt versus p Using Baseline Parameters')

    plt.axhline(y=0, color="k", linestyle="--", label="dpdt=0")
    plt.ylim((-5,10))
    plt.xlabel(f'p53 concentration (uM)')
    plt.ylabel("dp/dt (uM/min)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return dpdt_sol

if __name__ == "__main__":
    print("1D diffeq script was called explicitly")