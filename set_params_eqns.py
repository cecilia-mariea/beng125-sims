import numpy as np

# parameters

# params used in 1D report
PARAMETERS_1D = {
    "kp" : 0.8, # uM / min
    "kp_p53" : 6, # uM/ min
    "dp" : 0.1, # 1 / min
    "Kp" : 2, # uM
    "Km" : 1,  # uM
    "lam" : 3, # # 1 / min
    "n" : 4, # dimensionless , cooperative binding
    "m": 2 # dimensionless, cooperative binding
}

# change these params to change baseline params throughout scripts
EXP_PARAMETERS_2D = {
    "kp" : 0.8, # p53 basal synthesis
    "kp_p53" : 6, # autoregulation of p53
    "dp" : 0.1, # natural degredation of p53
    "Kp" : 2, # half saturation constant for TP53
    "Km" : 1,  # half saturation constant for MDM2 mRNA
    "lam" : 3, # ubiquitination of p53
    "n" : 4, #  cooperativity of TP53 binding
    "m": 2, # cooperativity of MDM2 mRNA binding
    "km": 0.5, # MDM2 basal synthesis
    "km_p53" : 1.5, # activation of MDM2
    "dm" : 0.2 # natural degredation of MDM2
}

# 1D eqns

def hill(p, n, K):
    p = np.asarray(p, dtype=float)
    return np.power(p, n) / (np.power(K, n) + np.power(p, n))

# 1D eqnsset_params_eqns.py
def dp(p, kp, kp_p53, dp, Kp, Km, lam, n, m):
    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * hill(p, m, Km)) * p)

def dpdt_1D(t, p, params):
    return dp(p, **params)

# 2D eqns

def dpdt_2D(vars, kp, kp_p53, dp, Kp, Km, lam, n, m, km, km_p53, dm):

    p, M = vars

    return (kp + kp_p53 * hill(p,n, Kp)) - ((dp + lam * M) * p)

def dMdt(vars, kp, kp_p53, dp, Kp, Km, lam, n, m, km, km_p53, dm):
    
    p, M = vars

    return (km + km_p53 * hill(p, m, Km)) - (dm * M)

def p53_MDM2(t, vars, params):

    return (dpdt_2D(vars, **params), dMdt(vars, **params))

if __name__ == "__main__":
    print("module loaded directly")