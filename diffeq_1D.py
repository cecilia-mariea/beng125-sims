import numpy as np
import matplotlib.pyplot as plt

BASE_PARAMETERS = {
    "kp" : 0.8, # uM / min
    "kp_p53" : 6, # uM/ min
    "dp" : 0.1, # 1 / min
    "Kp" : 2, # uM
    "Km" : 1,  # uM
    "lam" : 3, # # 1 / min
    "n" : 4, # dimensionless , cooperative binding
    "m": 2 # dimensionless, cooperative binding
} 

def dpdt(p, kp, kp_p53, dp, Kp, Km, lam, n, m):

    def hill(p, n, K):
        return (p ** n) / (K ** n + p ** n)

    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * hill(p, m, Km)) * p)

def dpdt_vs_p(dpdt, base_parameters, param_label="", param_vals=[]):

    p_vals = np.linspace(-2.5,5,200)
    plt.figure(figsize=(8,5))

    if bool(param_label):

        for val in param_vals: 

            params = base_parameters.copy()
            params[param_label] = val

            dpdt_sol = np.array([dpdt(p, **params) for p in p_vals])

            plt.plot(p_vals, dpdt_sol, label=f"{param_label}={val}", alpha=0.7)
            plt.title(f'dp/dt versus p, Varying {param_label}')

    else: 

        dpdt_sol = np.array([dpdt(p, **base_parameters) for p in p_vals] )

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

if __name__ == "__main__":

    # vary lambda

    dpdt_vs_p(dpdt, BASE_PARAMETERS, "lam", np.linspace(1.75, 3, 5))
    dpdt_vs_p(dpdt, BASE_PARAMETERS, "kp_p53", np.linspace(6, 14, 5))