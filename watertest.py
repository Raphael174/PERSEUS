import numpy as np
import matplotlib.pyplot as plt

## injector parameter

#test parameter
Dt_water_start = 10**(-7)
pas_Dt = 10**(-8)
diff_Re = 5 #percentage of error on the Re

#injector geometry
n = 4 #number of injector
D0 = 6*10**-3
Ds = 10*10**-3
Dt_inj = 10**(-3)
At_inj = n*np.pi*(Dt_inj**2)/4 #total area of tangential port
mass_flow_f = 0.61
mass_flow_o = 1.65

# water data
muWat=0.001
rhoWat=997
sigmaWat= 72*10**(-3)

# Lox data
rhoo=1141 #kg/m3 LOX
muo=451.36*10**(-6)  #Pa s LOX
sigmao=0.0132 #N/m

# ethanol data
rho = 789
muL=10**-3
sigmaL = 21*10**(-3)

delta_P = 3*10**5

## fonctions
def velocity(mass_flow,rho,At): # avec le Cd
    return (mass_flow)/(rho*At)

def Reynolds(rho,Dt,v,mu):
    return rho*Dt*v/mu

def Dt(Re,rho,mu,v):
    return mu*Re/(rho*v)

def CdLefebvre(At,D0,Ds):
    return 0.35*np.sqrt(At/(Ds*D0))*((Ds/D0)**0.25)

def water_test(Re_test):
# initialisation
    Dt_water = Dt_water_start
    At_water = n*np.pi*(Dt_water**2)/4
    mass_flow_w = np.sqrt((delta_P**2)*rhoWat*(CdLefebvre(At_water,D0,Ds)**2)*At_water)
    Re_w = Reynolds(rhoWat,Dt_water,velocity(mass_flow_w,rhoWat,At_water),muWat)

# search for the minimum Dt
    while np.abs(Re_w-Re_f) > diff_Re*Re_test/100:
        Dt_water += pas_Dt
        At_water = n*np.pi*(Dt_water**2)/4
        mass_flow_w = np.sqrt((delta_P**2)*rhoWat*(CdLefebvre(At_water,D0,Ds)**2)*At_water)
        Re_w = Reynolds(rhoWat,Dt_water,velocity(mass_flow_w,rhoWat,At_water),muWat)

    Dt_water_min = Dt_water
    mass_flow_w_min = np.sqrt((delta_P**2)*rhoWat*(CdLefebvre(At_water,D0,Ds)**2)*At_water)

#search for the maximum Dt
    while np.abs(Re_w-Re_f) < diff_Re*Re_f/100:
        Dt_water += pas_Dt
        At_water = n*np.pi*(Dt_water**2)/4
        mass_flow_w = np.sqrt((delta_P**2)*rhoWat*(CdLefebvre(At_water,D0,Ds)**2)*At_water)
        Re_w = Reynolds(rhoWat,Dt_water,velocity(mass_flow_w,rhoWat,At_water),muWat)

    Dt_water_max = Dt_water - pas_Dt
    At_water_max = n*np.pi*(Dt_water_max**2)/4
    mass_flow_w_max = np.sqrt((delta_P**2)*rhoWat*(CdLefebvre(At_water_max,D0,Ds)**2)*At_water_max)

#show result
    print('Dt : [',Dt_water_min,';',Dt_water_max,'] at',diff_Re,'% (in m)')
    print('mass flow : [',mass_flow_w_min,';',mass_flow_w_max,'] at',diff_Re,'%')



## calculs

# real injector Reynolds
Re_f = Reynolds(rho,Dt_inj,velocity(mass_flow_f,rho,At_inj),muL)
Re_o = Reynolds(rhoo,Dt_inj,velocity(mass_flow_o,rhoo,At_inj),muo)

# water test values
print('ethanol')
water_test(Re_f)
print('\nLOX')
water_test(Re_o)