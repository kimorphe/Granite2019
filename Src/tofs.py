#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

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
       print("f1,f2=",f1,f2)

       self.C=np.sum(Count,axis=1)
    def show(self,ax,fsz=16):
        ext=[self.time[0],self.time[-1],-self.ycod[0],-self.ycod[-1]]
        Vmax=np.sum(self.C[:])/100
        #ax.imshow(self.C, extent=ext,aspect="auto",origin="lower",cmap="jet",interpolation="bilinear",vmin=0,vmax=Vmax)
        im=ax.imshow(self.C, extent=ext,aspect="auto",origin="lower",cmap="jet",interpolation="bilinear",vmin=0,vmax=1500)
        ax_div=make_axes_locatable(ax) 
        cax=ax_div.append_axes("right",size="5%",pad="2.5%")
        cbar=colorbar(im,cax=cax,orientation="vertical");
        cbar.ax.tick_params(labelsize=fsz-2)

        ax.set_xlim(ext[0:2])
        ax.set_ylim(ext[2:4])
        ax.set_ylabel("y [mm]",fontsize=fsz)
        ax.set_xlabel("time [$\mu$ sec]",fontsize=fsz)
        ax.tick_params(labelsize=fsz)
    def time_stats(self,deg=1):
        indx=np.argmax(self.C,axis=1)
        self.tmax=self.time[indx]

        #jndx=np.argmax(self.C,axis=0)
        #self.ymax=self.ycod[jndx]
        #print("ymax=",self.ymax)
        
        Coef0=np.polyfit(self.tmax,self.ycod,deg)
        Coef1=np.polyder(Coef0)
        P0=np.poly1d(Coef0)
        P1=np.poly1d(Coef1)

        self.yfit=P0(self.tmax)
        self.cy=P1(self.tmax)   # phase velocity (based on max prob. TOF)

        rr=self.yfit-self.ycod
        res=np.mean(rr*rr)
        Lfit=np.mean(self.yfit*self.yfit)
        Ldat=np.mean(self.ycod*self.ycod)
        print(" RMS =",res)
        print(" RMS(normalized)=",res/np.sqrt(Lfit*Ldat)*100,"%")
        #print("Correlation=",np.mean(self.yfit*self.ycod)/np.sqrt(Lfit*Ldat))
        print(" <c>_tmin=",np.mean(self.cy))

        [T,Y]=np.meshgrid(self.time,self.ycod)
        nsum=np.sum(self.C,axis=1)
        tave=np.sum(self.C*T,axis=1)/nsum

        nsumy=np.sum(self.C,axis=0)
        yave=np.sum(self.C*Y,axis=0)/nsumy
        yvar=np.sum(self.C*Y*Y,axis=0)/nsumy;
        self.ysig=np.sqrt(yvar-yave*yave)
        self.yave=yave  

        tvar=np.sum(self.C*T*T,axis=1)/nsum;
        self.tsig=np.sqrt(tvar-tave*tave)
        self.tave=tave  

        coef0=np.polyfit(self.tave,self.ycod,deg)
        coef1=np.polyder(coef0)
        p0=np.poly1d(coef0)
        p1=np.poly1d(coef1)
        print(" <c>_tave=",np.mean(self.cy))
        self.vyfit=p0(self.tave)
        self.vy=p1(self.tave)   # phase velocity (based on mean TOF)


if  __name__=="__main__":

    H=Hist()    # Histogram data (y,w,t)
    S=Slice()   # Slice and stack (marginalize) Histogram over [w,w+dw]
    H.load("hist_ywt.dat")

    fsz=12  # label size

    DG=2    # degree of polynomial for curve fitting

    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    Stot=Slice()    
    Stot.set(H) # get the whole Histogram data
    Stot.show(ax) # show stacked histogram
    Stot.time_stats(deg=DG) # obtain statisics
    ycod=-Stot.ycod
    ax.plot(Stot.tmax,ycod,"w-")   # max. probability TOF curve (data)
    ax.plot(Stot.tmax,-Stot.yfit,"k-")   # max. probablitiy TOF curve (fitted) 
    ax.plot(Stot.tave,ycod,"y-")   # mean TOF curve (data)
    ax.plot(Stot.tave-Stot.tsig,ycod,"m--") # mean TOF+stdev
    ax.plot(Stot.tave+Stot.tsig,ycod,"m--") # mean TOF-stdev

    #ax.plot(Stot.time,-Stot.ymax,"oy")
    ax.plot(Stot.time,-Stot.yave,"ow")
    ax.plot(Stot.time,-Stot.yave+Stot.ysig,"--w")
    ax.plot(Stot.time,-Stot.yave-Stot.ysig,"--w")
    print("yvar=",Stot.ysig)


    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    S.set(H,fmin=0.8,fmax=1.2)  # get Hsitogram for given frequency band
    S.show(bx)  # show marginalized histogram
    S.time_stats() # obtain statistics concerning TOF 
    bx.plot(S.tmax,S.ycod,"w-") # TOF(data)
    bx.plot(S.tmax,S.yfit,"k-") # TOF(fitted)
    bx.plot(S.tave,S.ycod,"y-") # mean TOF
    bx.plot(S.tave-S.tsig,S.ycod,"m--") # mean TOF + stdev
    bx.plot(S.tave+S.tsig,S.ycod,"m--") # mean TOF - stdev

    fig3=plt.figure()
    fig4=plt.figure()
    fig5=plt.figure()
    cx1=fig3.add_subplot(111)
    cx2=fig4.add_subplot(111)
    cx3=fig5.add_subplot(111)
    Df=0.1  # frequncy band width [MHz] 
    Nf=int((H.freq[-1]-H.freq[0])/Df) # number of slices 
    Df=(H.freq[-1]-H.freq[0])/Nf    # frequency bandwidth (adjusted)
    fs=np.arange(Nf)*Df+H.freq[0]   # segmented frequencies
    cy=[]
    vy=[]
    for freq in fs:
        S.set(H,fmin=freq,fmax=freq+Df) # get Histogram data
        S.time_stats(deg=DG)    # obtain TOF statistics
        #txt=str(freq)+"MHz"
        #cx0.plot(-S.ycod,S.tave,label=txt)     # TOF(ave) as a function of ycod
        #cx0.plot(-S.ycod,S.tmax,"--",label=txt) # TOF(max prob.) as a function of ycod
        #cx3.plot(-S.ycod,S.tsig/S.tave) # normalized stdev{TOF} as a function of ycod
        cx3.plot(-S.ycod,S.tsig) # normalized stdev{TOF} as a function of ycod
        cy.append(-np.mean(S.cy))   # spatial average velocity (from tmax)
        vy.append(-np.mean(S.vy))   # spatial average velocity (from tave)
        cx2.plot(-S.vyfit,-S.vy,"k")
        cx2.plot(-S.yfit,-S.cy,"g")
    #cx3.plot(-Stot.ycod,Stot.tsig/S.tave,"k--",linewidth=2)
    cx3.plot(-Stot.ycod,Stot.tsig,"k--",linewidth=2)
    cx1.grid(True)
    cx2.grid(True)
    cx3.grid(True)
    fs=fs+0.5*Df
    cx1.plot(fs,cy,"-s",markersize=8,color="g",label="max prob") # plot velocity (max prob. TOF) as a function of freq.
    cx1.plot(fs,vy,"-v",markersize=8,color="k",label="mean") # plot velocity (ave TOF) as a function of freq. 
    cx1.set_ylim([2.5,3.5])
    cx2.set_ylim([2.0,4.0])

    cx1.set_xlabel("frequency [MHz]",fontsize=fsz)
    cx1.set_ylabel("phase velocity [km/sec]",fontsize=fsz)
    plt.legend()

    ax.tick_params(labelsize=fsz)
    bx.tick_params(labelsize=fsz)
    cx1.tick_params(labelsize=fsz)
    cx2.tick_params(labelsize=fsz)
    cx3.tick_params(labelsize=fsz)
    cyb=-np.mean(Stot.cy)
    vyb=-np.mean(Stot.vy)
    print("Total cy=",np.mean(Stot.cy))
    print("Total vy=",np.mean(Stot.vy))
    fmin=H.freq[0];
    fmax=H.freq[-1];
    cx1.hlines(cyb,fmin,fmax,colors="g",linestyles="dashed")
    cx1.hlines(vyb,fmin,fmax,colors="k",linestyles="dashed")
    cx1.set_xlim([fmin,fmax])

    cx2.set_xlim([-S.ycod[0],-S.ycod[-1]])
    cx2.set_xlabel("x [mm]",fontsize=fsz)
    cx2.set_ylabel("phase velocity [km/sec]",fontsize=fsz)

    cx3.set_xlim([-S.ycod[1],-S.ycod[-1]])
    cx3.set_xlabel("x [mm]",fontsize=fsz)
    cx3.set_ylabel("standard deviation (normalized) ",fontsize=fsz)
    plt.legend()

    fig6=plt.figure()
    cx4=fig6.add_subplot(111)
    cx4.plot(Stot.time,Stot.ysig,"o-")
    cx4.grid(True)

    plt.show()

    fig1.savefig("tof_all.png",bbox_inches="tight")
    fig2.savefig("tof_fbnd.png",bbox_inches="tight")
    fig3.savefig("vels_w.png",bbox_inches="tight")
    fig4.savefig("vels_y.png",bbox_inches="tight")
    fig5.savefig("stdev_TOF.png",bbox_inches="tight")


