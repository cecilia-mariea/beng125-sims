# imports
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# set base parameters, change if needed

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


# define ODE with auto regulation
def dpdt(t, p, kp, kp_p53, dp, Kp, Km, lam, n, m):

    def hill(p, n, K):
        return (p ** n) / (K ** n + p ** n)

    return (kp + kp_p53 * hill(p,n,Kp)) - ((dp + lam * hill(p, m, Km)) * p)

def solve_ODE_stable(t_span, t_eval, p0, params):
    return solve_ivp(
        fun = lambda  t, p: dpdt(t, p, **params), 
        t_span = t_span,
        y0 = [p0],
        t_eval = t_eval,
        rtol=1e-8, # relative tolerances for better asymptotic behavior
        atol=1e-10 # absolute tolerances for better asymptotic behavior
    )

def solve_ODE_unstable(t_span, t_eval, p0, params):
    return solve_ivp(
        fun = lambda t, p: -1 * dpdt(t, p, **params),
        t_span = t_span,
        y0 = [p0],
        t_eval = t_eval,
        rtol=1e-8, # relative tolerances for better asymptotic behavior
        atol=1e-10 # absolute tolerances for better asymptotic behavior
    )

# run sims and plot
def t_sim(t0, tf, p0, params, dt=1000):

    """
    
    One 1D simulations for p53 dynamics for different initial values of p53 concentration.

    t0 (int) : initial time for simulations
    tf (int) : final time for simulations
    p0 (tuple) : the initial conditions for simulations, in format (start, stop, step)  
    params (dict) : 
    dt (int) : number of time steps, defaulted at 1000

    """

    # set up simulation conditions
    t_span = (t0, tf)
    t_eval = np.linspace(t_span[0], t_span[1], dt)
    p0 = np.linspace(*p0)
    ss = []

    plt.figure(figsize=(8,5))
    for init_cond in p0: 
        stable_sol = solve_ODE_stable(
            t_span,
            t_eval,
            init_cond,
            params
        )
        ss.append(stable_sol.y[0,-1])
        plt.plot (stable_sol.t, stable_sol.y[0], label=f'p53(t=0)= {init_cond:.2f}')

    # plot
    plt.xlabel('Time (min)')
    plt.ylabel('Concentration of p53 (uM)')
    plt.title('Dynamics of p53 concentration for different initial conditions')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.tight_layout()
    plt.grid(True)
    plt.show()
    
    return ss

def main():
    # change simulation conditions
    t0 = 0
    tf = 5 
    min_p = 0.3 
    max_p = 1.2 
    num_init_cond = 10
    return t_sim(t0, tf, (min_p, max_p, num_init_cond), BASE_PARAMETERS, dt=500)

if __name__ == '__main__': 
    sol = main()
