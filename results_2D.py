import numpy as np

# 2D solutions
def figure6():

    from solve_2D import pairwise, t_sim

    # if many init cond, pass a iterable of init cond pairs (p,M)

    p0 = np.linspace(1,5,5)
    M0 = np.linspace(1,5,5)

    # init_cond = pairwise(p0, M0)
    t0 = 0
    tf = 10 

    sols = t_sim(t0, tf, (2,4))

# phase plane, example trajectory, nullclines
def figure7():

    from plot_phase_plane import plot_phase_plane
    from solve_2D import solve_sys

    p_range = (0, 4)
    M_range = (0, 4)    
    grid_size = 100 
    t0 = 0
    tf = 20 
    t_span = (t0, tf)
    t_eval = np.linspace(*t_span, 1000)
    init_cond = (4,4)

    sol = solve_sys(t_span, t_eval, init_cond)

    plot_phase_plane(p_range, M_range, grid_size, include_nullcline=True)

# param sweeping
def figure8():
   pass 