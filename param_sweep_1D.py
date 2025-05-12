from itertools import combinations
from math import ceil
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from bifurcations_1D import dpdt, find_fixed_points 

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

def bistability_analysis(dpdt, param_dict, base_params, init_cond):

    lables = list(param_dict.keys())
    values = list(param_dict.values())

    bistability_map = np.zeros((len(values[0]), len(values[1])))

    for i, x in enumerate(values[0]):
        for j, y in enumerate(values[1]):
            parameters = base_params.copy()
            parameters[lables[0]] = x 
            parameters[lables[1]] = y 

            fixed_points = find_fixed_points(dpdt, parameters, init_cond)

            bistability_map[i, j] = len(fixed_points)
                
    return bistability_map

def plot_all_bistability_data(bistability_data):

    num_plots = len(bistability_data)
    n = ceil(np.sqrt(num_plots))
    fig, axes = plt.subplots(n, n, figsize=(12,12))
    axes = axes.flatten()

    for ax in axes[num_plots:]:
        ax.axis('off')

    for (param_pair, data), ax in zip(bistability_data.items(), axes): 

        param1_vals = param_ranges[param_pair[0]]
        param2_vals = param_ranges[param_pair[1]]


        im = ax.imshow(data, 
                       origin="lower", 
                       aspect="auto", 
                       cmap="coolwarm",
                       extent=[param1_vals[0], param2_vals[-1], param1_vals[0], param1_vals[-1]],
                       vmin=1,
                       vmax=3
                       )

        fig.colorbar(im, ax=ax, label="# of Distinct Steady States")
        ax.set_title(f"{param_pair[0]} vs {param_pair[1]}")
        ax.set_xlabel(param_pair[0])
        ax.set_ylabel(param_pair[1])

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    param_ranges = {
        "kp" : np.linspace(0, 5, 50), 
        "kp_p53" : np.linspace(0, 5, 50),
        "dp" : np.linspace(0, 5, 50), 
        "Kp" : np.linspace(0.1, 5, 50),
        "Km" : np.linspace(0.1, 5, 50),  
        "lam" : np.linspace(0, 5, 50),
    }

    pairs = list(combinations([param for param in param_ranges], 2))
    bistability_data ={} 
    init_cond = np.linspace(0.1, 10, 20)

    for pair in tqdm(pairs):
        param_dict = {
            pair[0] : param_ranges[pair[0]],
            pair[1] : param_ranges[pair[1]]
        } 

        bistability_data[pair] = bistability_analysis(dpdt, param_dict, BASE_PARAMETERS, init_cond)

    plot_all_bistability_data(bistability_data)