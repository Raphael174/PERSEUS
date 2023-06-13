"""
Useful equations
"""
import numpy as np

def velocity2mdot(Vj, rho, do):
    #Vj : single oifice jet velocity
    #rho: liquid propellant density
    #do : orifice diameter

    d2A = np.pi * do**2 /4
    return Vj * rho * d2A

def mdot2velocity(mdot, rho, do):
    #Vj : single oifice jet velocity
    #rho: liquid propellant density
    #do : orifice diameter

    d2A = np.pi * do**2 /4
    return mdot/(rho * d2A)

def Reynolds(rho_L, mu_L, V, D):
    return rho_L*V*D/mu_L

def Webber(rho_L, Uj, do, sigma):
    return rho_L*Uj**2*do/sigma