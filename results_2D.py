import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

# 2D solutions
def figure6():

    from solve_2D import pairwise, t_sim

    # if many init cond, pass a iterable of init cond pairs (p,M)

    t0 = 0
    tf = 100 

    init_cond=np.array([(2,4), (4,2)])

    sols = t_sim(t0, tf, (1,1))

# phase plane, example trajectory, nullclines
def figure7():

    from plot_phase_plane import plot_phase_plane
    from solve_2D import solve_sys

    p_range = (0, 10)
    M_range = (0, 5)   
    grid_size = 100 
    t0 = 0
    tf = 20 
    t_span = (t0, tf)
    t_eval = np.linspace(*t_span, 1000)
    init_cond = (4,4)

    sol = solve_sys(t_span, t_eval, init_cond)

    plot_phase_plane(p_range, M_range, grid_size, include_nullcline=True)

# bifurcation digram

def figure8():

    from bifurcations_2D import calc_fixed_pts, classify_fixed_pt

    km_vals = np.linspace(0,14,50)
    kp_p53_vals = np.linspace(0,14,50)
    
    a_stable, b_stable = [], []
    a_unstable, b_unstable = [], []
    a_cycle, b_cycle = [], []

    for km in tqdm(km_vals):
        for kp in kp_p53_vals: 

            sols = calc_fixed_pts("Km", "kp_p53", km, kp)

            classification = [classify_fixed_pt(pt, "Km", "kp_p53", km, kp, spirals=True) for pt in sols]

            if any(cl == "stable" for cl in classification):
                a_stable.append(km)
                b_stable.append(kp)

            elif any(cl == "unstable" for cl in classification):

                a_unstable.append(km)
                b_unstable.append(kp)
            
            elif any(cl == "cycle" for cl in classification):

                a_cycle.append(km)
                b_cycle.append(kp)

    plt.figure(figsize=(8,8))

    plt.scatter(a_stable, b_stable, marker="x", color="red", label="stable spirals")

    plt.scatter(a_unstable, b_unstable, marker="o", color="blue", label="unstable sprials")

    plt.scatter(a_cycle, b_cycle, marker="*", color="black", label="limit cycles" )

    plt.xlabel("Parameter Km")
    plt.ylabel("Parameter kp_p53")
    plt.title("Behavior of Attractors in Parameter Space")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
