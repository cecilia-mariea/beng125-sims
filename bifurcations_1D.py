# imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.differentiate import derivative

from set_params_eqns import dp as dpdt, PARAMETERS_1D as params


# find fixed points
def find_fixed_points(params, init_cond):

    def f(p):
        return dpdt(p, **params)

    steady_states = []

    for cond in init_cond:
        sol = root(f,cond)
        if sol.success and sol.x[0] > 0:
            if not any(np.isclose(sol.x[0], ss, rtol=1e-4) for ss in steady_states):
                steady_states.append(sol.x[0]) 
    
    return sorted(steady_states)

# classify stability of the fixed points that were found
def classify_stability(params, ss):
    
    def f(p):
        return dpdt(p, **params)

    dfdp = derivative(f, ss).df

    return "stable" if dfdp < 0 else "unstable" 

# draw bifuctation diagrams
def bifurcation_diagram(param_label, param_vals, init_cond):
    
    bifurcation_data = []

    for param in param_vals:
        params_copy = params.copy()
        params_copy[param_label] = param
        
        fixed_points = find_fixed_points(params_copy, init_cond)

        for pt in fixed_points:
            stability = classify_stability(params_copy, pt)
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
    print("1D bifurcations module called directly")