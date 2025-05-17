# imports
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from set_params_eqns import dpdt_1D as dpdt
from set_params_eqns import PARAMETERS_1D

def solve_ODE(t_span, t_eval, p0, params):
    return solve_ivp(
        fun = lambda  t, p: dpdt(t, p, params), 
        t_span = t_span,
        y0 = [p0],
        t_eval = t_eval,
        rtol=1e-8, # relative tolerances for better asymptotic behavior
        atol=1e-10 # absolute tolerances for better asymptotic behavior
    )

# run sims and plot
def t_sim(t0, tf, p0, params, dt=1000):

    # set up simulation conditions
    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], dt)
    p0 = np.linspace(*p0)
    ss = []

    plt.figure(figsize=(8,5))
    for init_cond in p0: 
        sol = solve_ODE(
            t_span,
            t_eval,
            init_cond,
            params
        )
        ss.append(sol.y[0,-1])
        plt.plot (sol.t, sol.y[0], label=f'p53(t=0)= {init_cond:.2f}')

    # plot
    plt.xlabel('Time (min)')
    plt.ylabel('Concentration of p53 (uM)')
    plt.title('Dynamics of p53 concentration for different initial conditions')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.grid(True)
    plt.show()
    
    return ss

if __name__ == '__main__': 
    t0 = 0
    tf = 5 
    min_p = 0.3 
    max_p = 1.2 
    num_init_cond = 10

    ss = t_sim(t0, tf, (min_p, max_p, num_init_cond), PARAMETERS_1D)
