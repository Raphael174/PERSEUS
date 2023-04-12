import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter    
import numpy as np
from ExperimentalSSI import *
from SingleSwirlCorrel import *
from ErgolData import *
from CorrelParam import *

#/////////////////////// PLOTS ////////////////////////
param=["Do","Lo","Ds","At","Ao","AP","mf","SMD","nt","alpha"]
#x=LoWei
##---PLOT THE SMD------
def plotSMD(ExpSMD,abscissa,name):
    i=-1
    Q=["SMDLEV","SMDJAS","SMDLIU","SMDWAN","SMDWEI"]
    plt.figure(1)
    for P in Q:
        i=i+1
        if "Separated" in name:
            x=np.linspace(1,len(ExpSMD[i]),len(ExpSMD[i]))
        else:
            x=abscissa
        plt.plot(x,ExpSMD[i],label=P)
    plt.xlabel(name)
    plt.ylabel('Sauter Mean Diameter')
    #plt.xscale('log')
    plt.legend()
    plt.savefig('PLOTS/SMDplot')
    plt.show()


##---PLOT THE CD---------


def plotCD(ExpCD,abscissa,name):
    i=-1
    Q=["CDLEV","CDFU1","CDFU2","CDABRAM","CDGIFF","CDWEI","CDHONG","CDLIU","CDJONES","CDBENJ"]
    Q=["CDLEV","CDGIFF","CDWEI","CDHONG","CDJONES","CDBENJ"]
    plt.figure(2)
    for P in Q:
        i=i+1
        if "Separated" in name:
            x=np.linspace(1,len(ExpCD[i]),len(ExpCD[i]))
        else:
            x=abscissa
        plt.plot(x,ExpCD[i],label=P)
    plt.xlabel(name)
    plt.ylabel('Discharge Coefficient')
    #plt.yscale('log')
    plt.legend()
    plt.savefig('PLOTS/dropletOscillation')
    plt.show()



##------PLOT THE AP---------


def plotAP(ExpAP,abscissa,name):
    i=-1
    Q=["APLEV","APFU1","APFU2","APABRAM", "APGIFF","APWei","APHONG","APJONES","APBENJ"]
    Q=["APLEV","APGIFF","APWei"]
    plt.figure(3)
    for P in Q:
        i=i+1
        if "Separated" in name:
            x=np.linspace(1,len(ExpAP[i]),len(ExpAP[i]))
        else:
            x=abscissa
        plt.plot(x,ExpAP[i],label=P)
    plt.xlabel(name)
    plt.ylabel('Pressure Drop')
    #plt.yscale('log')
    plt.legend()
    plt.savefig('PLOTS/APplot')
    plt.show()

##--------PLOT SMD EVOLUTION WITH PARAMETER--------------
def plotSMDlev(name,fvalue,svalue,AP,lowlim,highlim):
    if (name=="Do"):
        Do=np.linspace(lowlim,highlim,20)
        Ds=fvalue
        Dt=svalue 
        aP=AP
    if (name=="Ds"):
        Ds=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Dt=svalue
        aP=AP
    if (name=="Dt"):
        Dt=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Ds=svalue
        aP=AP
    if (name=="aP"):
        aP=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Ds=svalue
        Dt=AP
    
    Dx=np.linspace(lowlim,highlim,20) 
    smd=SMDlevRED(Do,Ds,Dt,aP)
    plt.xlabel(name)
    plt.ylabel('SMD')
    plt.plot(Dx,smd)
    plt.savefig('PLOTS/SMDplot'+ name)
    plt.show()

##--------PLOT SMD EVOLUTION WITH IMPINGING--------------
def plotSMDimp(name,fvalue,svalue,AP,lowlim,highlim,SMDimp,SMDimp2):
    if (name=="Do"):
        Do=np.linspace(lowlim,highlim,20)
        Ds=fvalue
        Dt=svalue 
        aP=AP
    if (name=="Ds"):
        Ds=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Dt=svalue
        aP=AP
    if (name=="Dt"):
        Dt=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Ds=svalue
        aP=AP
    if (name=="aP"):
        aP=np.linspace(lowlim,highlim,20)
        Do=fvalue
        Ds=svalue
        Dt=AP
    
    Dx=np.linspace(lowlim,highlim,20)
    smdImp=np.ones(len(Dx))*SMDimp
    smdImp2=np.ones(len(Dx))*SMDimp2
    smd=SMDlevRED(Do,Ds,Dt,aP)
    plt.xlabel(name)
    plt.ylabel('SMD')
    plt.plot(Dx,smd,Dx,smdImp,Dx,smdImp2)
    plt.savefig('PLOTS/SMDplot'+ name)
    plt.show()
