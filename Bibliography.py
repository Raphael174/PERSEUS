import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter    
import numpy as np
import random
from sympy import symbols, Eq, solve


from SingleSwirlCorrel import *   #Imports all the correlationns known for parameters [Basic, Cd,thickness,Vt,AP,Beta,BuL,SMD,mf_ev]
from SSAnimation import *         #Imports the animation file
from constants.ErgolData import *           #Imports the data from the Ergols [density, surface tension, diffussion,...]
from ExperimentalSSI import *     #Imports the experimental data from a set of injectors tested by "Mr.Wei" in his paper
from CorrelParam import *         #Imports the functions to estimate performance parameters with the best working correlations
from Plots import *               #Imports plot funcions to simplify the plotting of the correlations and experimental data.
# from BibliographyImpinging import*


#////////////////////APROXIMATIONS///////////////////////////

#Parameters that need to be optimised to get the best performance out of the injector.
#Input of the system

Do = 0.002;
Lo=0.02 #orifice length, the whole length of the air coloumn
Ds = 0.015;
Dt=0.004
nt=3
At = np.pi*(Dt/2)**2*nt
Ao= np.pi*(Do/2)**2
aP=30*10**5
#mf = 0.02;
alpha = np.pi/4 #convergence angle of the injector measured in the wall from inside, zero for an open end swirl
theta=0 #Angle of the inlet ports

#From the literature (Improved Semiempirical Correlation to Predict Sauter Mean Diameter for Pressure-Swirl Atomizers) by Xiao Wei and Huang Yong
#we have found several injectors that have been tested, and whose SMD have been measured.
#Function "getExpData" outputs a matrix with rows ["Do","Lo","Ds","At","Ao","AP","mf","SMD","nt","alpha"], (Not the same order as above)
#being the parameters of each injector and 18 columns corresponding to all the injectors tested.
#The function allows to order the data following a certain parameter.
#For example, getExpData("Do") would order all the injectors in increasing outlet diameter.

ExpData=getExpData("Do")    #Choose a string between the parameters apporximations to order the data accordingly
                            #In this case, data is ordered following an increasing Do



#////////////////////CORRELATED PARAMETERES////////////////////

#From the input parameters, using correlations defined in file "SingleSwirlCorrel", the rest of the parameters affecting the
#injector are estimated. Only one correlation is used for every parameter (The one that seems to work the best with the 
#experimental results). This way all the parameters are estimated in the file "CorrelParam". 


#Param=EstimParamAnimate(At,Ds,Do,mf,Ao,Dt,alpha,Lo) #This function creates the animation parameters, given the approximations
                                                    #that the user chooses 
                                                    #and the output is "[Do,HELEV,BETA,FRECUNS,LBU,SMD,MDROP,VT]"". Which is
                                                    #the input of the "animate" function defined in file "SSAnimation".

ExpSMD=ExpParmPlotSMD(ExpData)          #This functions create a matrix with rows ["SMDLEV","SMDJAS","SMDLIU","SMDWAN","SMDWEI"]
ExpCD=ExpParmPlotCD(ExpData)            #and columns the 18 injectors with the input values from the experiments of "Mr.Wei"
ExpAP=ExpParmPlotAP(ExpData)            #This way the correlations can be comparated with the real values to see which work better
                                        #As a conclusion, Levfebre's correlations tend to work pretty well

DO=getExpAbscissa("Do")           #function getExpAbscissa("strg") gets the abscissa needed for plotting. It can be any value
                                        #from ["Do","Lo","Ds","At","Ao","AP","mf","SMD","nt","alpha"]. If not any particular order is 
                                        #requested, "original" would output the order decided by "Mr.Wei" in his paper.

#////////////////////ANIMATION/////////////////////////

#animateSSI(Param)                       #Function from file "SSAnimation" that creates an injector simulation given the correlated
                                         #parameters obtained in the previous section. A screen capture of the animation is saved in 
                                         #folder PLOTS

#/////////////////// EXPERIMENTAL DATA PLOTS /////////////////

plotSMD(ExpSMD,DO,"Do") #Plot contains the ordinate, abscissa and name of the abscissa. If the name contains the word "Separate"
plotCD(ExpCD,DO,"Do")   #the plot will separate the abscissa if the same x value has different y values. The plots will be displayed
plotAP(ExpAP,DO,"Do")   #on screen and saved in folder PLOTS. As the same plot is done, it overrides in the folder. 

#From this plots we conclude that Levfebre's correlation is the best for all applications and that the range of validity is
#Do[0.001, 0.005]
#Ds[0.012  0.018]
#Dt[0.003  0.0045]

#We use the validated ranges to plot the corrleations
plotSMDlev("Ds",Do,Dt,aP,0.012,0.018)
plotSMDlev("Do",Ds,Dt,aP,0.001,0.005)
plotSMDimp("Do",Ds,Dt,aP,0.001,0.0156,SMDImpinging,SMDImpinging2)
plotSMDlev("Dt",Do,Ds,aP,0.003,0.0045)
plotSMDlev("aP",Do,Ds,Dt,10*10**5,40*10**5)

#Do must be minimized and is the most relevant parameter
#Ds must be maximized and its variation does not affect much (Follow the literature)
#Dt must be minimized but does not have such an impacta as Do
#AP must be maximized but fixed by tank pressure
print(SMDlevRED(0.001,0.018,0.003,aP))                                  #Target of 29.8 microns - Target 50 microns following NASA - 100 microns theory


