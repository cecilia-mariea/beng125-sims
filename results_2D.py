import numpy as np

# figure 6 2D solutions
def figure6():
    from plot_phase_plane import plot_phase_plane
    from solve_2D import solve_sys

    p_range = (0, 10)
    M_range = (0, 10)    
    grid_size = 10 

    t0 = 0
    tf = 15
    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    sol = solve_sys(t_span, t_eval, (0,1))

    plot_phase_plane(p_range, M_range, grid_size, traj=sol, include_traj=True, normalized=True)

# phase plane

# nullclines

# param sweeping