# imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.cm import ScalarMappable

# custom imports
from set_params_eqns import  dMdt 
from set_params_eqns import dpdt_2D as dpdt, PARAMETERS_2D as params

# calculate trajectories
# calculate 2D vector field
def plot_phase_plane(p_range, M_range, grid_size, traj="", include_traj=False, normalized=False):

    # set up points in the vector field
    p = np.linspace(p_range[0], p_range[1], grid_size)
    m = np.linspace(M_range[0], M_range[1], grid_size)
    P, M = np.meshgrid(p,m)

    # set up arrays to populate with specific solutions
    dpdt_spec_sol = np.zeros_like(P)
    dMdt_spec_sol = np.zeros_like(M)

    # solve for the specific solution at all points in the vector feild
    for i in range(grid_size):
        for j in range(grid_size):

            vars = (P[i,j], M[i,j])

            dpdt_spec_sol[i,j] = dpdt(vars, **params)
            dMdt_spec_sol[i,j] = dMdt(vars, **params)
    
    # find the cartesian magnitude
    magnitudes = np.sqrt(dpdt_spec_sol**2 + dMdt_spec_sol**2)

    if normalized:
        dpdt_spec_sol = dpdt_spec_sol / magnitudes
        dMdt_spec_sol = dMdt_spec_sol / magnitudes
    

    plt.figure(figsize=(8,8))

    plt.quiver(P,M, dpdt_spec_sol, dMdt_spec_sol, angles='xy', scale_units='xy', scale=1, color='b')

    plt.xlim(p_range[0], p_range[1])
    plt.ylim(M_range[0], M_range[1])
    plt.xlabel("p53 concentration (uM)")
    plt.ylabel("MDM2 concentration (uM)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.title("P53 and MDM2 Phase Plane")

    if include_traj:
        start_point = round(traj.y[0][0],2), round(traj.y[1][0],)
        plt.plot(traj.y[0], traj.y[1], 'r-', label=f"t=0, {start_point}")
        plt.plot(start_point[0], start_point[1], 'ro')
        
    plt.legend()
    plt.show()

    return P, M, dpdt_spec_sol, dMdt_spec_sol, magnitudes

if __name__ == '__main__':
    print("2D phase plane script was explicitly called")