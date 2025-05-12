# imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.misc import derivative

# parameters

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

# find fixed points
def find_fixed_points(dpdt, params, init_cond):

    steady_states = []

    def f(p):
        return dpdt(p, **params)

    for cond in init_cond:
        sol = root(f,cond)
        if sol.success and sol.x[0] > 0:
            if not any(np.isclose(sol.x[0], ss, rtol=1e-4) for ss in steady_states):
                steady_states.append(sol.x[0]) 
    
    return sorted(steady_states)

# classify stability of the fixed points that were found
def classify_stability(dpdt, ss, params):
    
    def f(p):
        return dpdt(p, **params)

    dfdp = derivative(f, ss, dx=1e-6)

    return "stable" if dfdp < 0 else "unstable" 

# draw bifuctation diagrams
def bifurcation_diagram(dpdt, param_label, param_vals, params, init_cond):
    
    bifurcation_data = []

    for param in param_vals:
        params_copy = params.copy()
        params_copy[param_label] = param
        
        fixed_points = find_fixed_points(dpdt, params_copy, init_cond)

        for pt in fixed_points:
            stability = classify_stability(dpdt, pt, params_copy)
            bifurcation_data.append((param, pt, stability))

    # seperate unstable and stable branches
    stable_param = [data[0] for data in bifurcation_data if data[2] == "stable"]
    unstable_param = [data[0] for data in bifurcation_data if data[2] == "unstable"]

    stable_pt = [data[1] for data in bifurcation_data if data[2] == "stable"]
    unstable_pt = [data[1] for data in bifurcation_data if data[2] == "unstable"]

    # plot 
    plt.figure(figsize=(8,5))
    plt.scatter(stable_param, stable_pt, color="blue", s=20, label="Stable Branch")

    plt.scatter(unstable_param, unstable_pt, color="red", s=20, label="Unstable Branch")

    plt.xlabel(f'{param_label}')
    plt.ylabel('p53 concentration (uM)')
    plt.title(f'Bifurcation Diagram: Fixed points versus {param_label}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    return bifurcation_data

if __name__ == "__main__":

    init_cond = np.linspace(1, 10, 20)

    # vary lambda 

    lam_vals = np.linspace(1,3.5, 100) 
    bifurcation_diagram(dpdt, 'lam', lam_vals, BASE_PARAMETERS, init_cond)

    # vary kp_p53

    kp_p53_vals = np.linspace(5,18, 100)
    bifurcation_diagram(dpdt, 'kp_p53', kp_p53_vals, BASE_PARAMETERS, init_cond)
