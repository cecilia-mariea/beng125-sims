# imports
import numpy as np
import matplotlib.pyplot as plt

# custom imports
from set_params_eqns import  dMdt, calc_nullcline
from set_params_eqns import dpdt_2D as dpdt, PARAMETERS_2D as params
from fixed_points_2D import calc_fixed_pts, classify_fixed_pt

# calculate 2D vector field
def plot_phase_plane(p_range, M_range, grid_size, normalized=True, include_nullcline=False):
    """
    `traj` should be a solve_ivp object if `include_traj` is True, can use solve_2D `solve_sys`
    """

    # set up points in the vector field
    p = np.linspace(p_range[0], p_range[1], grid_size)
    m = np.linspace(M_range[0], M_range[1], grid_size)
    P, M = np.meshgrid(p,m)

    # set up arrays to populate with specific solutions
    dpdt_spec_sol = np.zeros_like(P)
    dMdt_spec_sol = np.zeros_like(M)

    # solve for the specific solution at all points in the vector field
    for i in range(grid_size):
        for j in range(grid_size):

            vars = (P[i,j], M[i,j])

            dpdt_spec_sol[i,j] = dpdt(vars, **params)
            dMdt_spec_sol[i,j] = dMdt(vars, **params)
    
    if normalized:
        # find the cartesian magnitude
        magnitudes = np.sqrt(dpdt_spec_sol**2 + dMdt_spec_sol**2)

        dpdt_spec_sol = dpdt_spec_sol / magnitudes
        dMdt_spec_sol = dMdt_spec_sol / magnitudes

    plt.figure(figsize=(6,6))

    plt.streamplot(P,M, dpdt_spec_sol, dMdt_spec_sol, color="orange", density=0.8)

    plt.xlim(p_range[0], p_range[1])
    plt.ylim(M_range[0], M_range[1])
    plt.xlabel("p (uM)")
    plt.ylabel("M (uM)")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.title("P53 and MDM2 Phase Plane")

    if include_nullcline:
        
        P, p_nullcline_M, M_nullcline_M = calc_nullcline(p, **params)

        plt.plot(P, p_nullcline_M, label="P53 nullcline", color="r", linewidth=2)
        plt.plot(P, M_nullcline_M, label="MDM2 nullcline", color="b", linewidth=2)

        fixed_pt = calc_fixed_pts()
        print(fixed_pt)

        classification = [classify_fixed_pt(pt) for pt in fixed_pt]
        print(classification)

        for pt in fixed_pt:
            x,y = pt
            plt.scatter(x,y,marker="o", color="black",s=40,label="fixed point")

    plt.legend(loc='upper left')
    plt.show()

    return P, M, dpdt_spec_sol, dMdt_spec_sol

if __name__ == '__main__':
    print("2D phase plane script was explicitly called")