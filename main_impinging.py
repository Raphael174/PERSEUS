"""
Testing correlations and plottingg graphs
scientific plotting: https://ranocha.de/blog/colors/

#   TODO: scientific plots
    TODO: unlike impinging injectors doublet
    TODO: optimal FOF OFO plot correlations, optiaml designs
    TODO: impact of 3D printing roughness
"""
#%%

from constants.ErgolData import *
from Impinging.doublets.D_avg_semi_empirical_correlations import D_drop_mean_LinStab, D_drop_mean_John1995
from Impinging.doublets.D_avg_correlations import D_drop_mean_Sweeney2016, D_drop_mean_Yoon2004, D_drop_mean_Lourme1986
from Impinging.useful.pressure_loss import DeltaPpipe
from Impinging.useful.useful_equations import mdot2velocity
import numpy as np

import matplotlib.pyplot as plt

#%%

OF = 1.4
mdot_tot = 2.04
mdot_ox = OF/(1+OF) * mdot_tot
print("mdot_tot_ox = ", mdot_ox)
N_likedoublets_ox = 2
varrying_mass_flow_lox = np.linspace(0.7*mdot_ox, 1.3*mdot_ox, 10)
#2 orifices in 1 pair * 6 ox elements
unitary_lox = varrying_mass_flow_lox/(2*N_likedoublets_ox)

#take d_o 3mm
d_o = 2e-3
Lo___do = 3
L_o = Lo___do * d_o
#assume
Lj = d_o*4

beta = 90 #deg

unitary_mdot_to_V = mdot2velocity(unitary_lox, rhoo, d_o)
Us_momentumBalance = unitary_mdot_to_V * np.cos(beta/2 * np.pi/180)

#%%
print("Range of orifice velocities")
print(unitary_mdot_to_V)
print("Resulting sheet velocities")
print(Us_momentumBalance)


""" COMPUTING DROP SIZE """
#%% Linear instability semi-empirical correlation 
# Initiated by the work of Dombrowski and expanded by John

mean_drop_size_LinStab = []
mean_drop_size_LinStab = D_drop_mean_LinStab(d_o, beta, Us_momentumBalance, sigmao, rhog, rhoo)*1e6

print("Estimated mean drop diameter from Linear Instability µm")
print(mean_drop_size_LinStab)

#Correlations
drop_size_john1995 = D_drop_mean_John1995(d_o, beta, Us_momentumBalance, sigmao, rhog, rhoo)*1e6
drop_size_micron_Sweeney2016 = D_drop_mean_Sweeney2016(d_o, rhoo, sigmao, Us_momentumBalance, beta)
#Yoon assumes bbeta = 90°
drop_size_Yoon2004 = D_drop_mean_Yoon2004(d_o, rhoo, rhog, sigmao, Us_momentumBalance)*1e6
drop_size_Lourme1986 = D_drop_mean_Lourme1986(d_o, rhog, muo, sigmao, Us_momentumBalance,Lj)
print(drop_size_Lourme1986)
# %%Pressure drop from pipe friction loss
Delta_P = []
Delta_P = DeltaPpipe(rhoo, muo, Us_momentumBalance, L_o, d_o)/1e5
print("Pressure loss bar")
print(Delta_P)



# %% plot results

normalize_over_orifice = d_o*1e6
#plt.scatter(Delta_P, mean_drop_size_LinStab, label = " Linear instability")
plt.scatter(Delta_P, drop_size_john1995, label = "John1995")
plt.scatter(Delta_P, drop_size_micron_Sweeney2016, label="Sweeney2016" )
plt.scatter(Delta_P, drop_size_Yoon2004, label="Yoon2004")
plt.scatter(Delta_P, drop_size_Lourme1986, label="Lourme1986")
plt.xlabel(r'$\Delta$ P [bar]')
plt.ylabel('Mean Drop Diameter [µm]')
plt.title("Mean drop size computations")
plt.grid()
plt.legend()
#plt.savefig("Like_Doublet_impinging_CorrelationComparison.png", dpi=400)
plt.show()

#%%