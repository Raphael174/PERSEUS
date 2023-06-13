import numpy as np

def D_drop_mean_Sweeney2016(d_o, rho_L, sigma, Us, theta):
    We = rho_L*Us**2*d_o/sigma

    return 4751*We**-0.38 * np.sin(theta*np.pi/180)**-0.46

def D_drop_mean_Yoon2004(d_o, rho_L, rho_g, sigma, Us):
    We = rho_L*Us**2*d_o/sigma
    rho_g__rho_L = rho_g/rho_L

    return d_o * 1.64 * rho_g__rho_L**(-1/6)*We**(-1/3)

def D_drop_mean_Lourme1986(d_o, rho_g, mu_L, sigma, Us, Lj):
    termV = (Us/30)**-0.95
    termD = (d_o/2)**0.3
    termL = (Lj/5)**-0.08
    termRhog= (rho_g/5)**-0.2
    termSigma = (sigma/(7.35*10**-2))**0.5
    termMu = (mu_L*1e3)**-0.04

    return 150*termV*termD*termL*termRhog*termSigma*termMu