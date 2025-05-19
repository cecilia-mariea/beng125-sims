# ubiquitous imports :|
import numpy as np

# figure 2
def figure2():

    from diffeq_1D import dpdt_vs_p

    p_vals = np.linspace(-2.5,5,200)
    dpdt_vs_p(p_vals) 


# figure 3
def figure3():

    from solve_1D import t_sim

    t0 = 0
    tf = 5 
    min_p = 0.3 
    max_p = 1.2 
    num_init_cond = 10

    ss = t_sim(t0, tf, (min_p, max_p, num_init_cond))

# figure 4
def figure4():

    from diffeq_1D import dpdt_vs_p 

    p_vals = np.linspace(-2.5,5,200)

    dpdt_vs_p(p_vals,"lam", np.linspace(1.75, 3, 5))
    dpdt_vs_p(p_vals,"kp_p53", np.linspace(6, 14, 5))

# figure 5
def figure5():

    from bifurcations_1D import bifurcation_diagram

    init_cond = np.linspace(1, 10, 20)

    # vary lambda 

    lam_vals = np.linspace(1,3.5, 100) 
    bifurcation_diagram('lam', lam_vals, init_cond)

    # vary kp_p53

    kp_p53_vals = np.linspace(5,18, 100)
    bifurcation_diagram('kp_p53', kp_p53_vals,init_cond)