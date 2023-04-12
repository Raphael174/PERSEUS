#////////////////////WEI EXPERIMENTAL DATA//////////////////////////
import numpy as np
DoWei=np.array([0.001, 0.002, 0.003, 0.001, 0.002, 0.003, 0.001, 0.002, 0.003,0.001,0.002,0.003,0.0004,0.0004,0.0004,0.0006,0.0006,0.0006])
LoWei=np.array([1,0.5,0.33,1,0.5,0.33,1,0.5,0.33,1,0.5,0.33,0.5,1.25,2.5,5,6.67,8.33])
DsWei=np.array([0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.018,0.012,0.012,0.012,0.002,0.002,0.002,0.002,0.002,0.002])
AtWei=np.array([8,8,8,12,12,12,16,16,16,12,12,12,0.59,0.59,0.59,0.59,0.59,0.59])*10**(-6)
AoWei=np.array([0.785,3.14,7.07,0.785,3.14,7.07,0.785,3.14,7.07,0.785,3.14,7.07,0.126,0.126,0.126,0.126,0.126,0.126])
APWei=np.ones(len(DoWei))*0.6*10**6
mfWei=np.array([0.013,0.032,0.0580,0.014,0.035,0.063,0.015,0.042,0.072,0.018,0.043,0.063,0.003,0.0029,0.0034,0.0059,0.00492,0.00477])
SMDWei=np.array([108,105,119,98,106,140,124,130,158,88,116,151,45,49,54,114,124,161])
ntWei=np.ones(len(DoWei))*4
alphaWei=np.ones(len(DoWei))*np.pi/4
ExpData=np.zeros((10,len(DoWei)))
ExpData[0]=DoWei
ExpData[1]=LoWei
ExpData[2]=DsWei
ExpData[3]=AtWei
ExpData[4]=AoWei
ExpData[5]=APWei
ExpData[6]=mfWei
ExpData[7]=SMDWei
ExpData[8]=ntWei
ExpData[9]=alphaWei

def getExpData(string):
    if(string=="Do"):
        arr=DoWei
    if(string=="Lo"):
        arr=LoWei
    if(string=="Ds"):
        arr=DsWei
    if(string=="At"):
        arr=AtWei        
    if(string=="Ao"):
        arr=AoWei
    if(string=="AP"):
        arr=APWei
    if(string=="mf"):
        arr=mfWei
    if(string=="SMD"):
        arr=SMDWei
    if(string=="nt"):
        arr=ntWei
    if(string=="alpha"):
        arr=alphaWei
    if(string=="original"):
        return ExpData

    arr2=np.array(np.sort(arr))
    ConnectM=np.zeros((len(arr),len(arr2)))
    Repeated=np.zeros(len(arr))
    saved=0;
    rep=0;
    alr=0;
    ii=-1
    jj=-1
    for i in arr2:
        ii=ii+1
        jj=-1
        alr=0
        if(i==saved):
            rep=rep+1
        else:
            rep=0
            saved=i
        for j in arr:
            jj=jj+1
            if (j==i) and (alr==rep):
                ConnectM[ii,jj]=1
            else:
                ConnectM[ii,jj]=0
            if(j==i):
                alr=alr+1
            
    ConnectM=np.linalg.inv(ConnectM)
    #print(np.dot(arr,ConnectM)-arr2)
    output=np.dot(ExpData,ConnectM)
    #print(output[0])
    #print(output[1])
    return output

def getExpAbscissa(string):
    if(string=="Do"):
        arr=DoWei
    if(string=="Lo"):
        arr=LoWei
    if(string=="Ds"):
        arr=DsWei
    if(string=="At"):
        arr=AtWei        
    if(string=="Ao"):
        arr=AoWei
    if(string=="AP"):
        arr=APWei
    if(string=="mf"):
        arr=mfWei
    if(string=="SMD"):
        arr=SMDWei
    if(string=="nt"):
        arr=ntWei
    if(string=="alpha"):
        arr=alphaWei
    if(string=="original"):
        return np.linspace(1,len(DoWei),len(DoWei))
    arr2=np.array(np.sort(arr))
    return arr2

"""""
class Injector:
    def __init__(self,Do,Lo,Ds,At,Ao,AP,mf,SMD) -> None:
        self.Do=Do
        self.Lo=Lo
        self.Ds=Ds
        self.At=At
        self.AP=AP
        self.Ao=Ao
        self.mf=mf
        self.SMD=SMD

    def getDo(self):
        return self.Do
    
    def getLo(self):
        return self.Lo
    
    def getDs(self):
        return self.Ds
    
    def getAt(self):
        return self.At
    
    def getAP(self):
        return self.AP
    
    def getAo(self):
        return self.Ao
    
    def getmf(self):
        return self.mf
    
    def getSMD(self):
        return self.SMD
    

    def InjSort(self, Injector):
        if Injector.Lo/Injector.Do>self.Lo/self.Do:
            return 1
        if Injector.Lo/Injector.Do<self.Lo/self.Do:
            return -1
        if Injector.Lo/Injector.Do==self.Lo/self.Do:
            return 0

Inj0=Injector(0.001,1,0.018,8,0.785,0.6*10**6,0.013,108)
Inj1=Injector(0.002,0.5,0.018,8,3.14,0.6*10**6,0.032,105)
Inj2=Injector(0.003,0.33,0.018,8,7.07,0.6*10**6,0.058,119)
Inj3=Injector(0.001,1,0.018,12,0.785,0.6*10**6,0.014,98)
Inj4=Injector(0.002,0.5,0.018,12,3.14,0.6*10**6,0.035,106)
Inj5=Injector(0.003,0.33,0.018,12,7.07,0.6*10**6,0.063,140)
Inj6=Injector(0.001,1,0.018,16,0.785,0.6*10**6,0.015,124)
Inj7=Injector(0.002,0.5,0.018,16,3.14,0.6*10**6,0.042,130)
Inj8=Injector(0.003,0.33,0.018,16,7.07,0.6*10**6,0.072,158)
Inj9=Injector(0.001,1,0.012,12,0.785,0.6*10**6,0.018,88)
Inj10=Injector(0.002,0.5,0.012,12,3.14,0.6*10**6,0.043,116)
Inj11=Injector(0.003,0.33,0.012,12,7.07,0.6*10**6,0.063,151)
Inj12=Injector(0.0004,0.5,0.002,0.59,0.126,0.6*10**6,0.003,45)
Inj13=Injector(0.0004,1.25,0.002,0.59,0.126,0.6*10**6,0.0029,49)
Inj14=Injector(0.0004,2.5,0.002,0.59,0.126,0.6*10**6,0.0034,54)
Inj15=Injector(0.0006,5,0.002,0.59,0.126,0.6*10**6,0.0059,114)
Inj16=Injector(0.0006,6.67,0.002,0.59,0.126,0.6*10**6,0.00492,124)
Inj17=Injector(0.0006,8.33,0.002,0.59,0.126,0.6*10**6,0.00477,161)


Injectors=[Inj0,Inj1,Inj2,Inj3,Inj4,Inj5,Inj6,Inj7,Inj8,Inj9,Inj10,Inj11,Inj12,Inj13,Inj14,Inj15,Inj16,Inj17]


"""""
