import numpy as np

# 2D solutions
def figure6():

    from solve_2D import pairwise, t_sim

    # if many init cond, pass a iterable of init cond pairs (p,M)

    p0 = np.linspace(1,5,5)
    M0 = np.linspace(1,5,5)

    init_cond = pairwise(p0, M0)
    t0 = 0
    tf = 15

    sols = t_sim(t0, tf, init_cond, init_range=True)

# phase plane and example trajectory
def figure7():
    from plot_phase_plane import plot_phase_plane
    from solve_2D import solve_sys

    p_range = (0, 5)
    M_range = (0, 5)    
    grid_size = 10 

    t0 = 0
    tf = 20 
    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    init_cond = (4,4)
    sol = solve_sys(t_span, t_eval, init_cond)

    plot_phase_plane(p_range, M_range, grid_size, sol, include_traj=True)

# nullclines

# param sweeping