#LOSSES in pipe

#%%
import numpy as np


#%%
#take
roughness = 0.0025e-3 #mm

def Altshul_II_friction(Re, do, roughness = roughness):
    r = roughness
    return (1.8*np.log10(Re/(0.135*Re*(r/do) + 6.5)))**(-2)


def DeltaPpipe(rho_L, mu_L, Vj, Lo, do):
    Re = rho_L*Vj*do/mu_L #Reynolds(rho_L, mu_L, Vj, do)
    f = Altshul_II_friction(Re, do)
    return 0.5 * rho_L * f * Vj**2 * Lo/do

# %%