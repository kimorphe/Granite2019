import numpy as np
import matplotlib.pyplot as plt

class Hist:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        dat=fp.readline()
        yrng=list(map(float,dat.strip().split(",")))
        print("yrng=",yrng)

        fp.readline()
        dat=fp.readline()
        frng=list(map(float,dat.strip().split(",")))
        print("frng=",frng)

        fp.readline()
        dat=fp.readline()
        trng=list(map(float,dat.strip().split(",")))
        print("trng=",trng)

        fp.readline()
        dat=fp.readline()
        Nd=list(map(int,dat.strip().split(",")))
        print("Nd=",Nd)

        fp.readline()
        amp=[]
        for row in fp:
            amp.append(float(row))
        amp=np.array(amp)
        print(np.shape(amp))
        amp=np.reshape(amp,Nd)

        self.ycod=np.linspace(yrng[0],yrng[1],Nd[0])
        self.freq=np.linspace(frng[0],frng[1],Nd[1])
        self.time=np.linspace(trng[0],trng[1],Nd[2])
        self.Ny=Nd[0];
        self.Nf=Nd[1];
        self.Nt=Nd[0];
        self.Count=amp;
class Slice:
    def set(self,H,fmin=-1,fmax=-1):
       self.ycod=H.ycod 
       self.time=H.time 
       self.Ny=H.Ny
       self.Nt=H.Nt
       if fmin==-1:
           fmin=H.freq[0];
       if fmax==-1:
           fmax=H.freq[-1];
       f1=fmin
       f2=fmax;

       nf1=np.argmin(np.abs(H.freq-f1));
       nf2=np.argmin(np.abs(H.freq-f2));
       f1=H.freq[nf1]
       f2=H.freq[nf2]
       self.freq=H.freq[nf1:nf2+1]
       Count=H.Count[:,nf1:nf2+1,:]
       print("f1=",f1)
       print("f2=",f2)
       print(np.shape(Count))

       self.C=np.sum(Count,axis=1)
    def show(self,ax):
        ext=[self.time[0],self.time[-1],self.ycod[0],self.ycod[-1]]
        Vmax=np.sum(self.C[:])/100
        #ax.imshow(self.C, extent=ext,aspect="auto",origin="lower",cmap="jet",interpolation="bilinear",vmin=0,vmax=Vmax)
        ax.imshow(self.C, extent=ext,aspect="auto",origin="lower",cmap="jet",interpolation="bilinear",vmin=0,vmax=1500)
    def time_stats(self):
        indx=np.argmax(self.C,axis=1)
        self.tmax=self.time[indx]
        
        deg=2
        Coef0=np.polyfit(self.tmax,self.ycod,deg)
        Coef1=np.polyder(Coef0)
        P0=np.poly1d(Coef0)
        P1=np.poly1d(Coef1)

        self.yfit=P0(self.tmax)
        self.cy=P1(self.tmax)
        print("c(y)=",self.cy)
        print("<c>=",np.mean(self.cy))

        rr=self.yfit-self.ycod
        res=np.sum(rr*rr)

        Lfit=np.sum(self.yfit*self.yfit)
        Ldat=np.sum(self.ycod*self.ycod)
        print("res(err[%])=",res,res/np.sqrt(Lfit*Ldat)*100,np.sum(self.yfit*self.ycod)/np.sqrt(Lfit*Ldat))

        [T,Y]=np.meshgrid(self.time,self.ycod)
        print("shape T",np.shape(T))
        print("shape Y",np.shape(Y))
        
        nsum=np.sum(self.C,axis=1)
        tave=np.sum(self.C*T,axis=1)/nsum
        print(tave)

        tvar=np.sum(self.C*T*T,axis=1)/nsum;
        self.tsig=np.sqrt(tvar-tave*tave)
        print("sig=",self.tsig)
        self.tave=tave




if  __name__=="__main__":
    H=Hist()
    S=Slice()
    H.load("hist_ywt.dat")

    S.set(H,fmin=1.8,fmax=2.0)

    D=np.sum(H.Count,axis=1)
    indx=np.argmax(D,axis=1)
    print(np.shape(D))
    print(np.shape(indx))


    fig=plt.figure()
    ax=fig.add_subplot(111)
    print("All frequency")
    Stot=Slice()
    Stot.set(H)
    Stot.show(ax)
    Stot.time_stats()
    ax.plot(Stot.tmax,Stot.ycod,"w-")
    ax.plot(Stot.tmax,Stot.yfit,"k-")
    ax.plot(Stot.tave,Stot.ycod,"y-")
    ax.plot(Stot.tave-Stot.tsig,Stot.ycod,"m--")
    ax.plot(Stot.tave+Stot.tsig,Stot.ycod,"m--")

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    S.show(bx)
    print("f in [0.6,1.0]")
    S.time_stats()
    bx.plot(S.tmax,S.ycod,"w-")
    bx.plot(S.tmax,S.yfit,"k-")
    bx.plot(S.tave,S.ycod,"y-")
    bx.plot(S.tave-S.tsig,S.ycod,"m--")
    bx.plot(S.tave+S.tsig,S.ycod,"m--")


    plt.show()
