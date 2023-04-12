import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter    
import numpy as np
import random
import sqlalchemy
import mongoengine


#-----------------------------Basic Parameters-----------------------------------   
def Acte(Ds,Do,Dt,nt):
    return (Ds-Dt)/Do/Dt**2/nt

def Kcte(At,Ds,Dt,Do):
    return At/Do/(Ds-Dt)

def Re(rho,Vt,Dt,muL):
    return rho*Vt*Dt/muL

def We(rho,Vt,Dt,sigmaL):
    return rho*Vt**2*Dt/sigmaL

def gamma(Ds,Dt,Do):
    return (Ds-Dt)/Do

def Da(Do,Re,alpha,Ds,Dt,Lo): #Must have a range of validity as for alpha close to zero must yield a value
    return Do*0.338*(1-np.exp(-1.45*10**(-4)*Re))*alpha**0.073*(Do/Ds)**0.424*(Dt/Ds)**(-0.732)*(Lo/Ds)**(-0.252)

def Xratio(Da,Do):
    return Da/Do #Soltani


#---------------------Discharge Coefficients-----------------------------
def CdFu1(Acte): #OpenEnd
    return 0.4354/Acte**0.877

def CdFu2(At,Do,gamma): #OpenEnd
    return 0.19*(At/Do**2)**0.65*gamma**(-2.13)

def CdAbram(Acte):
    return 0.432/Acte**0.64

def CdGiff(Xratio):
    return 1.17*np.sqrt((1-Xratio)**3/(1+Xratio))

def CdLev(At,Ds,Do): #Ds -> Swirl diameter //// Do -> Output diameter //// At -> Total area for the input swirl channels
    return 0.35*np.sqrt(At/Ds/Do)*(Ds/Do)**0.25;

def CdHong(At,Do,Ds,Dt):
    return 0.44*(At/Do**2)**(0.84*((Ds-Dt)/Do)**(-0.52))*((Ds-Dt)/Do)**(-0.59)

def CdLiu(Acte,Lo,Do,Ds):
    return 0.721*(1/Acte)**0.416*(Lo/Do)**(-0.0558)*(Ds/Do)**0.147

def CdJones(Do,rhoL,Vo,muL,Lo,Ls,Ds,At):
    return 0.45*(Do*Vo*rhoL/muL)**(-0.02)*(Lo/Do)**(-0.03)*(Ls/Ds)**0.05*(At/Ds/Do)**0.52*(Ds/Do)**0.23

def CdBenj(Do,rhoL,Vo,muL,Lo,Ls,Ds,At):
    return 0.466*(Do*Vo*rhoL/muL)**(-0.027)*(Lo/Do)**(-0.229)*(Ls/Ds)**0.091*(At/Ds/Do)**0.517*(Ds/Do)**0.187

def CdWei(AP,mf,Ao,rhoL):
    return mf/Ao/np.sqrt(2*rhoL*AP)

#--------------------------Out Film thickeness---------------------------
def helev(Do,mf,muL,rho,AP): #Convergend end theoretical
    return 3.66*(Do*mf*muL/rho/AP)**0.25

def heSuy(Do,mf,muL,rho,AP): #Open end exp
    return 2.7*(Do*mf*muL/rho/AP)**0.25

def heFu(Do,mf,muL,rho,AP): #Used in Double Coaxiaul Raph correlation Open end exp
    return 3.1*(Do*mf*muL/rho/AP)**0.25

def heKim(Do,mf,muL,rho,AP,Lc): #Open end exp
    return 3.1*(Do*mf*muL/rho/AP)**0.25*(Lc/Do)**0.6



#----------------------------- EXIT VELOCITY & AP ---------------------------------------
def AP(mf, Cd, At, rhoL):
    return mf**2/2/rhoL/(Cd*At)**2

def mf(AP,Cd,At,rhoL):
    return np.sqrt(AP*2*rhoL*Cd**2*At**2)

def Vtfun(mf,At,Do,he,rho):
    return mf/(At-np.pi*(Do/2-he)**2)/rho

#--------------------Spray Cone Angle--------------------------
def BetaLev(Xratio): #Theoretical
    return np.arccos(((1-Xratio)/(1+Xratio))**0.5)

def BetaGiff(CdGiff,Xratio): #Theoretical
    return np.arcsin(np.pi*CdGiff/2/(1+np.sqrt(Xratio))/np.sqrt(np.pi**2/32/Xratio**2*(1-Xratio)**3))

def BetaLiu(theta,Acte,Do,Lo,Ds): #Experimental Converge End
    return np.arccos(0.302*(1+np.tan(theta))**0.414*(1/Acte)**0.35*(Lo/Do)**0.043*(Ds/Do)**0.026+0.612)

def BetaFu(Acte,Re): #Experimental open end
    return np.arctan(0.033*Acte**0.338*Re**0.249)


#-----------------------BreakUp Length-----------------------
def LbuMoon(rhog,rho,Do,sigmaL,Vt,Beta,he,kunst): #Thinning Film
    return 1/2/np.sin(Beta)*((3/2*12*np.sqrt(2)*np.sin(Beta)*Vt*np.cos(Beta)*np.sqrt(he*(Do-he)/kunst/(rhog/rho*Vt**2-kunst*sigmaL/rho))))**(2/3)+(Do-he)**(3/2)

def LbuFu(rho,sigmaL,he,Beta,rhog,Vt): #Conical Film
    return 0.82*(rho*sigmaL*2.5*he*np.cos(Beta)/rhog**2/Vt**2)**0.5

def LbuIna(rho,sigmaL,he,Beta,rhog,Vt): #Conical Film
    return 0.2175*(rho*sigmaL*2.5*he*np.cos(Beta)/rhog**2/Vt**2)**0.3


#----------------------------SMD-------------------------------
def SMDlev (sigma,mu,mf,rhog,AP):
    return 2.25*sigma**0.25*mu**0.25*mf**(0.25)*rhog**(-0.25)*AP**(-0.5)*10**6

def SMDrad (sigma,muL,mf,rho,AP):
    return 7.3*sigma**0.6*muL**(0.2)*mf**(0.25)*rho**(-0.2)*AP**(-0.4)*10**6

def SMDjas (sigma,muL,mf,rho,AP):
    return 4.4*sigma**0.6*muL**(0.16)*mf**(0.22)*rho**(-0.16)*AP**(-0.43)*10**6

#def SMDball (At,muL,Do,AP):
#    return 0.436*muL**(0.55)*Do**(-0.05)*At**(-0.24)*AP**(-0.74)*10**6

#def SMDban (sigmaL,muL,mf,rhog,AP):
#    return (sigmaL*muL*mf/rhog)**0.25*3.29*10000*AP**(-0.5)

def SMDLiu (sigmaL,muL,rho,Do,rhog,AP,Ds,Lo,Beta):
    return 0.5536*sigmaL**0.25*muL**0.25*rho**0.125*Do**0.5*rhog**(-0.25)*AP**(-0.375)*10**6*(Ds/Do)**0.33*(Lo/Do)**0.122*(1+np.tan(Beta/2))    

def SMDWan (sigmaL,muL,rhog,AP,he,Beta,rho):
    return 4.52*(sigmaL*muL**2/rhog/AP**2)**0.25*(he*np.cos(Beta))**0.25+0.39*(sigmaL*rho/rhog/AP)**0.25*(he*np.cos(Beta))**0.75*10**6

def Dav (sigmaL,muL,rhog,AP,he,Beta,rho):
    return 2.11*(sigmaL*muL**2/rhog/AP**2)**0.25*(he*np.cos(Beta))**0.25+0.62*(sigmaL*rho/rhog/AP)**0.1*(he*np.cos(Beta))**0.9

#-----------------------------Droplet Evaporation----------------------
def mfdrop(rdrop,Diff,rho,Yko):
    return 4*np.pi*rdrop*rho*Diff*np.log(1/(1-Yko))

def mdrop(rdrop,rho):
    return 4/3*np.pi*rdrop**3*rho

def dist(mdrop,mfdrop,Vflight):
    return Vflight*mdrop/mfdrop


#----------------------OSCILLATING FLOW AT EXIT---------------------------------
def kunst(rhog,Vt,sigmaL): #Unstable wave number
    return rhog*Vt**2/2/sigmaL

def frecUns(rhog,Vt,sigmaL,rho,kunst,helev):
    return np.sqrt((2*rhog/rho*Vt**2*kunst**2-sigmaL*kunst**3/rho)/kunst/helev)

def h(Do,Beta,he,x):
    return he*(Do-he)/(Do-he+2*x*np.sin(Beta))

def Lbu(rhog,rho,Do,sigmaL,Vt,Beta,he,kunst):
    return 1/2/np.sin(Beta)*((3/2*12*np.sqrt(2)*np.sin(Beta)*Vt*np.cos(Beta)*np.sqrt(he*(Do-he)/kunst/(rhog/rho*Vt**2-kunst*sigmaL/rho))))**(2/3)+(Do-he)**(3/2)

def amplitude(amp0,x,frecUns,Vt):
    return amp0*np.exp(frecUns*x/Vt)

def oscillation(frecUns,t,x,amp0,Vt):
    return np.sin(frecUns*x+t)*amplitude(amp0,x,frecUns,Vt)


