import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter    
import numpy as np
import random
from SingleSwirlCorrel import *


class Star:    

    def __init__(self,x,y) -> None:
        self.x=x
        self.y=y
    
    def __init__(self,x,Beta,MAX,MDROP):
        self.mdrop=MDROP*10**12
        randnum=np.random.rand()*MAX*2-MAX
        self.x=x*np.cos(Beta)+np.sin(Beta)*(randnum)
        self.y=np.cos(Beta)*(randnum)-x*np.sin(Beta)
        self.mean=-x*np.sin(Beta)
    
    def move(self,t,Beta):
        self.x=self.x+t/5/10**3*np.cos(Beta)+t/5/10**3*50*np.sin(Beta)*(self.y-self.mean)/0.5
        self.y=self.y-t/5/10**3*np.sin(Beta)+t/5/10**3*50*np.cos(Beta)*(self.y-self.mean)/0.5
        self.mean=self.mean-t/5/10**3*np.sin(Beta)
        self.mdrop=np.abs(self.mdrop-self.mdrop/20)


    def getx(self):
        return self.x
    
    def gety(self):
        return self.y

    def getmdrop(self):
        return self.mdrop



#----------------------------------ANIMATION---------------------------

def animateSSI(Param):
    Do =        Param[0]
    HELEV=      Param[1]
    Beta=       Param[2]
    FRECUNS=    Param[3]
    LBU=        Param[4]
    MDROP=      Param[6]
    Vt=         Param[7]

    sizelist=[]
    SMDlist=[]
    xlist=[]
    ylist=[]

    discretH=10
    AMP0=h(Do,Beta,HELEV,0)*0.01


    def hAnim(x):
        return h(Do,Beta,HELEV,x)

    def oscilAnim(x,t):
        return np.sin(FRECUNS*x+t)*amplitude(AMP0,x,FRECUNS,Vt)


    fig, axs= plt.subplots(1,1)
    counter2=0
    for t in np.linspace(0,60,200):
        plt.clf()
        counter2=counter2+1
        for x in np.linspace(0,LBU,120):

            for i in np.linspace(0,hAnim(x),discretH):
                xlist.append(x*np.cos(Beta)+(oscilAnim(x,t)-hAnim(x)/2+i)*np.sin(Beta))
                ylist.append(-x*np.sin(Beta) + (oscilAnim(x,t)-hAnim(x)/2+i)*np.cos(Beta))
                sizelist.append(1)

        
        for i in SMDlist: i.move(120/200,Beta)
        for i in np.linspace(0,4,5):
            SMDlist.append(Star(LBU,Beta,amplitude(AMP0,LBU,FRECUNS,Vt),MDROP))
        for i in SMDlist: 
            xlist.append(i.getx())
            ylist.append(i.gety())
            #sizelist.append(np.abs((LBU-(i.getx()-LBU))*10**3))
            sizelist.append(i.getmdrop())




        plt.xlabel('downstream length')
        plt.ylabel('raidial length')
        plt.xlim(-0.002,0.02)
        plt.ylim(-0.015,0.002)
        plt.scatter(xlist,ylist,s=sizelist)
        plt.pause(0.001)
        if (counter2==190):
            plt.savefig('dropletEvolution')
        xlist=[]
        ylist=[]
        sizelist=[] 
