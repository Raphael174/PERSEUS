import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter    
import numpy as np
import random


from SingleSwirlCorrel import *
from SSAnimation import *
from ErgolData import *
from ExperimentalSSI import *

        
def EstimParamAnimate(At,Ds,Do,mf,Ao,Dt,alpha,Lo):
    CD=CdLev(At,Ds,Do)
    aP=AP(mf,CD,Ao,rho)
    HELEV=helev(Do,mf,muL,rho,aP)
    Vt=Vtfun(mf,At,Do,HELEV,rho)
    KUNST=kunst(rhog,Vt,sigmaL)
    FRECUNS=frecUns(rhog,Vt,sigmaL,rho,KUNST,HELEV)

    RE=Re(rho,Vt,Dt,muL)
    DA=Da(Do,RE,alpha,Ds,Dt,Lo)
    XRATIO=Xratio(DA,Do)
    BETA=BetaLev(XRATIO)
    LBU=Lbu(rhog,rho,Do,sigmaL,Vt,BETA,HELEV,KUNST)

    SMD=SMDlev(sigmaL,muL,mf,rho,aP)
    MFDROP = mfdrop(SMD*10**(-6)/2,DiffKer,rhog,Yko)
    MDROP = mdrop(SMD*10**(-6)/2,rho)
    VT= Vtfun(mf,At,Do,HELEV,rho)

    return [Do,HELEV,BETA,FRECUNS,LBU,SMD,MDROP,VT]

def ExpParmPlotSMD(ExpData):
    DoWei=ExpData[0]
    LoWei=ExpData[1]
    DsWei=ExpData[2]
    AtWei=ExpData[3]
    AoWei=ExpData[4]
    APWei=ExpData[5]
    mfWei=ExpData[6]
    SMDWei=ExpData[7]
    ntWei=ExpData[8]
    alphaWei=ExpData[9]
    
    DtWei=np.sqrt(AtWei/np.pi)
    heWei=helev(DoWei,mfWei,muWat,rhoWat,APWei)
    VtWei=Vtfun(mfWei,AtWei,DoWei,heWei,rhoWat)
    ReWei=Re(rhoWat,VtWei,DtWei,muWat)
    DaWei=Da(DoWei,ReWei,alphaWei,DsWei,DtWei,LoWei)
    XratioWei=Xratio(DaWei,DoWei)
    BETAWei=BetaLev(XratioWei)
    AoWei=np.pi*DoWei**2/4

    SMDLEV=SMDlev(sigmaWat,muWat,mfWei,rhoWatg,APWei)
    SMDJAS=SMDjas(sigmaWat,muWat,mfWei,rhoWat,APWei)
    SMDLIU=SMDLiu(sigmaWat,muWat,rhoWat,DoWei,rhoWatg,APWei,DsWei,LoWei,BETAWei)
    SMDWAN=SMDWan(sigmaWat,muWat,rhoWatg,APWei,heWei,BETAWei,rhoWat)
    result=np.array( [SMDLEV,SMDJAS,SMDLIU,SMDWAN,SMDWei])
    return result

def ExpParmPlotCD(ExpData):
    DoWei=ExpData[0]
    LoWei=ExpData[1]
    DsWei=ExpData[2]
    AtWei=ExpData[3]
    AoWei=ExpData[4]
    APWei=ExpData[5]
    mfWei=ExpData[6]
    SMDWei=ExpData[7]
    ntWei=ExpData[8]
    alphaWei=ExpData[9]

    DtWei=np.sqrt(AtWei/np.pi)
    heWei=helev(DoWei,mfWei,muWat,rhoWat,APWei)
    VtWei=Vtfun(mfWei,AtWei,DoWei,heWei,rhoWat)
    ReWei=Re(rhoWat,VtWei,DtWei,muWat)
    DaWei=Da(DoWei,ReWei,alphaWei,DsWei,DtWei,LoWei)
    XratioWei=Xratio(DaWei,DoWei)
    BETAWei=BetaLev(XratioWei)
    AoWei=np.pi*DoWei**2/4

    CDLEV=CdLev(AtWei,DsWei,DoWei)
    ActeWei=Acte(DsWei,DoWei,DtWei,ntWei)
    CDFU1=CdFu1(ActeWei)
    gammaWei=gamma(DsWei,DtWei,DoWei)
    CDFU2=CdFu2(AtWei,DoWei,gammaWei)
    CDABRAM=CdAbram(ActeWei)
    CDGIFF=CdGiff(XratioWei)
    CDWEI=CdWei(APWei,mfWei,AoWei,rhoWat)
    CDHONG=CdHong(AtWei,DoWei,DsWei,DtWei)
    CDLIU=CdLiu(ActeWei,LoWei,DoWei,DsWei)
    CDJONES=CdJones(DoWei,rhoWat,VtWei,muWat,LoWei,LoWei/2,DsWei,AtWei)
    CDBENJ=CdBenj(DoWei,rhoWat,VtWei,muWat,LoWei,LoWei/2,DsWei,AtWei)
    result=np.array( [CDLEV, CDGIFF,CDWEI,CDHONG,CDJONES,CDBENJ])
    return result


def ExpParmPlotAP(ExpData):
    DoWei=ExpData[0]
    LoWei=ExpData[1]
    DsWei=ExpData[2]
    AtWei=ExpData[3]
    AoWei=ExpData[4]
    APWei=ExpData[5]
    mfWei=ExpData[6]
    SMDWei=ExpData[7]
    ntWei=ExpData[8]
    alphaWei=ExpData[9]

    DtWei=np.sqrt(AtWei/np.pi)
    heWei=helev(DoWei,mfWei,muWat,rhoWat,APWei)
    VtWei=Vtfun(mfWei,AtWei,DoWei,heWei,rhoWat)
    ReWei=Re(rhoWat,VtWei,DtWei,muWat)
    DaWei=Da(DoWei,ReWei,alphaWei,DsWei,DtWei,LoWei)
    XratioWei=Xratio(DaWei,DoWei)
    BETAWei=BetaLev(XratioWei)
    AoWei=np.pi*DoWei**2/4

    CDLEV=CdLev(AtWei,DsWei,DoWei)
    ActeWei=Acte(DsWei,DoWei,DtWei,ntWei)
    CDFU1=CdFu1(ActeWei)
    gammaWei=gamma(DsWei,DtWei,DoWei)
    CDFU2=CdFu2(AtWei,DoWei,gammaWei)
    CDABRAM=CdAbram(ActeWei)
    CDGIFF=CdGiff(XratioWei)
    CDWEI=CdWei(APWei,mfWei,AoWei,rhoWat)
    CDHONG=CdHong(AtWei,DoWei,DsWei,DtWei)
    CDLIU=CdLiu(ActeWei,LoWei,DoWei,DsWei)
    CDJONES=CdJones(DoWei,rhoWat,VtWei,muWat,LoWei,LoWei/2,DsWei,AtWei)
    CDBENJ=CdBenj(DoWei,rhoWat,VtWei,muWat,LoWei,LoWei/2,DsWei,AtWei)

    APLEV=AP(mfWei,CDLEV,AoWei,rhoWat)
    APFU1=AP(mfWei,CDFU1,AoWei,rhoWat)
    APFU2=AP(mfWei,CDFU2,AoWei,rhoWat)
    APABRAM=AP(mfWei,CDABRAM,AoWei,rhoWat)
    APGIFF=AP(mfWei,CDGIFF,AoWei,rhoWat)
    APHONG=AP(mfWei,CDHONG,AoWei,rhoWat)
    APBENJ=AP(mfWei,CDBENJ,AoWei,rhoWat)
    APJONES=AP(mfWei,CDJONES,AoWei,rhoWat)
    result=np.array( [APLEV,APGIFF,APWei])
    return result

def SMDlevRED(Do,Ds,Dt,aP):
    At=np.pi*Dt**2/4
    Cd=CdLev(At,Ds,Do)
    Ao1=np.pi*Do**2/4
    mf1=mf(aP,Cd,Ao1,rho)
    return SMDlev(sigmaL,muL,mf1,rhog,aP)