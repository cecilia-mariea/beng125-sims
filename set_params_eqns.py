import numpy as np

""" 
these are most of the things that are hard-coded in for the entire module, all functionality should be preserved if changes are made to this file
"""

# parameters

# params used in 1D report
PARAMETERS_1D = {
    "kp" : 0.8, # p53 basal synthesis
    "kp_p53" : 6, # autoregulation of p53
    "dp" : 0.1, # natural degredation of p53
    "Kp" : 2, # half saturation constant for TP53
    "Km" : 1,  # half saturation constant for MDM2 mRNA
    "lam" : 3, # # ubiquitination of p53
    "n" : 4, # cooperativity of TP53 binding
    "m": 2 # cooperativity of MDM2 mRNA binding
}

# params used in 2D report 
PARAMETERS_2D = {
    "kp": 0.1,         # weak basal p53 synthesis
    "kp_p53": 6,     # strong positive feedback on p53
    "dp": 0.1,          # moderate p53 degradation
    "Kp": 2.0,          # p53 half-saturation for self-activation
    "Km": 1,          # Km for MDM2 induction, use Km=6 for unstable spiral
    "lam": 3.0,         # strong MDM2-mediated p53 degradation
    "n": 4,             # sharp p53 feedback (high cooperativity)
    "m": 3,             # sharp MDM2 activation
    "km": 0.01,         # weak basal MDM2 synthesis
    "km_p53": 3.0,      # strong MDM2 induction by p53
    "dm": 0.3           # relatively slow MDM2 degradation
}

# 1D eqns

def hill(p, n, K):
    p = np.asarray(p, dtype=float)
    return np.power(p, n) / (np.power(K, n) + np.power(p, n))

# 1D
def dp(p, kp, kp_p53, dp, Kp, Km, lam, n, m):
    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * hill(p, m, Km)) * p)

dpdt_1D = lambda t, p, params: dp(p, **params)

# 2D eqns

def dpdt_2D(vars, kp, kp_p53, dp, Kp, Km, lam, n, m, km, km_p53, dm):

    p, M = vars

    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * M) * p)

def dMdt(vars, kp, kp_p53, dp, Kp, Km, lam, n, m, km, km_p53, dm):
    
    p, M = vars

    return (km + km_p53 * hill(p, m, Km)) - (dm * M)

p53_MDM2_dt = lambda t, vars, params: (dpdt_2D(vars, **params), dMdt(vars, **params))

p53_MDM2 = lambda vars, params: (dpdt_2D(vars,**params), dMdt(vars, **params))

# nullclines : solve for M when dpdt = 0 etc. and hard code in

def calc_nullcline(p, kp, kp_p53, dp, Kp, Km, lam, n, m, km, km_p53, dm):

    p_nullcline_M = np.zeros_like(p)
    M_nullcline_M = np.zeros_like(p)

    def p_nullcline(p):

        if p > 0: 
            numerator = (kp + kp_p53 * hill(p, n, Kp)) - dp * p
            return numerator / (lam * p)
        else: 
            return np.nan

    def M_nullcline(p):
        return (km + km_p53 * hill(p, m, Km)) / dm

    for i, p_val in enumerate(p):

        p_nullcline_M[i] = p_nullcline(p_val) 
        M_nullcline_M[i] = M_nullcline(p_val)

    return p, p_nullcline_M, M_nullcline_M

if __name__ == "__main__":
    print("set params and eqns script loaded directly")