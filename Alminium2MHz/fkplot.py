#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

class BNDL:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        dat=fp.readline()
        Xa=list(map(float,dat.strip().split(",")))
        print("Xa=",Xa)

        fp.readline()
        dat=fp.readline()
        dx=list(map(float,dat.strip().split(",")))
        print("dx=",dx)

        fp.readline()
        dat=fp.readline()
        Nd=list(map(int,dat.strip().split(",")))
        print("Nd=",Nd)

        fp.readline()
        dat=fp.readline()
        dat=list(map(float,dat.strip().split(",")))
        f1=dat[0]
        df=dat[1]
        Nf=int(dat[2])
        print("f1,df,Nf=",f1,df,Nf)
        #Nf=int(Nf/2)

        fp.readline()
        amp=[]
        for row in fp:
            dat=row.strip().split(",")
            z=float(dat[0])+1j*float(dat[1])
            amp.append(z)
        amp=np.array(amp)
        Nx=Nd[0]
        Ny=Nd[1]
        #amp=amp.astype(np.float)
        print(np.size(amp))
        print(Nx*Ny*Nf)
        amp=np.reshape(amp,[Nx,Ny,Nf])
        print(np.shape(amp))
        self.amp=amp
        self.Nx=Nx
        self.Ny=Ny
        self.Nf=Nf
        self.dky=1/dx[1]/Ny;
        self.ky=np.arange(Ny)*self.dky
        self.freq=np.arange(Nf)*df+f1;
        self.df=df;
        self.Xa=Xa
        self.dx=dx
        self.xcod=np.arange(Nx)*dx[0]+Xa[0]
        self.ycod=np.arange(Ny)*dx[1]+Xa[1]
        fp.close()
    def get_index(self,val,axis):
        if axis==0:
            indx=np.argmin(np.abs(self.xcod-val))
        if axis==1:
            indx=np.argmin(np.abs(self.ycod-val))
        if axis==2:
            indx=np.argmin(np.abs(self.freq-val))
        if axis==3:
            indx=np.argmin(np.abs(self.ky-val))
        return(indx)
    def get_cod(self,indx,axis):
        cod=0.0;
        if axis==0:
            indx=indx%self.Nx;
            cod=self.xcod[indx]
        if axis==1:
            indx=indx%self.Ny;
            cod=self.ycod[indx]
        if axis==2:
            indx=indx%self.Nf;
            cod=self.freq[indx]
        if axis==2:
            indx=indx%self.Ny;
            cod=self.ky[indx]
        return(cod)
    def get_domain(self,axis):
        x=-self.xcod
        y=-self.ycod
        f= self.freq
        ky=-self.ky
        if axis==0:
            ext=[f[0],f[-1],y[0],y[-1]]
        if axis==1:
            ext=[f[0],f[-1],y[0],y[-1]]
        if axis==2:
            ext=[y[0],y[-1],x[0],x[-1]]
        if axis==3: # for fk-plot
            ext=[f[0],f[-1],ky[0],ky[-1]]
        return(ext)
    def stack(self,ax):
        S=np.mean(self.amp,axis=ax)
        return(S)



if __name__=="__main__":

    dir_name="./Scopes"
    fname="scopes_win.fft"
    fname=dir_name+"/"+fname

    bndl=BNDL()
    bndl.load(fname)

    fsz=12
    fig1=plt.figure(figsize=(6,4))
    ax=fig1.add_subplot(111)
    # Mean B-scan
    S=bndl.stack(0)
    FK=np.fft.fft(S,axis=0)
    ext=bndl.get_domain(3)
    findx=np.argmax(np.abs(FK),axis=1)
    f_peak=bndl.freq[findx] # peak frequency
    kindx=np.argmax(np.abs(FK),axis=0)
    k_peak=bndl.ky[kindx]   # peak wave number 
    FK=FK/np.max(np.abs(FK[:])) # normalize f-k spectrum

    ky_max=-0.6
    f_max=2.5
    imax=bndl.get_index(ky_max,3)
    jmax=bndl.get_index(f_max,2)

    #Df=0.15
    Df=0.06
    inc=int(Df/bndl.df)

    im=ax.imshow(np.abs(FK),extent=ext,interpolation="none",aspect="auto",cmap="jet",origin="lower")
    indx=np.arange(0,jmax+1,inc)
    #ax.plot(bndl.freq[indx],-k_peak[indx],".y",markersize="6")
    axdiv=make_axes_locatable(ax)
    cax=axdiv.append_axes("right",size="7%",pad="2%")
    cb=colorbar(im,cax=cax)


    fs=bndl.freq[0:jmax+1];

    c=-bndl.freq[indx]/k_peak[indx]

    deg=1
    coef=np.polyfit(k_peak[0:jmax+1],fs,deg)
    pk=np.poly1d(coef)
    pkd=np.poly1d(np.polyder(coef))
    #ax.plot(pk(k_peak[0:jmax+1]),-k_peak[0:jmax+1],"w-")
    cg1=-pkd(k_peak[indx]);
    print("k-f Linfit coef=",coef)

    kfit=k_peak[indx]
    ffit=pk(k_peak[indx])
    c_fit=ffit/kfit
    c_fit_mean=np.mean(c_fit[3:])
    print("c_fit=",c_fit)
    print("c_fit_mean=",c_fit_mean)

    deg=2
    coef=np.polyfit(k_peak[0:jmax+1],fs,deg)
    pk=np.poly1d(coef)
    pkd=np.poly1d(np.polyder(coef))
    #ax.plot(pk(k_peak[0:jmax+1]),-k_peak[0:jmax+1],"g--")
    cg2=-pkd(k_peak[indx]);

    ax.tick_params(labelsize=fsz)
    ax.set_xlim([0,3.0])
    ax.set_ylim([0,1.2])
    ax.set_xlabel("frequency [MHz]",fontsize=fsz);
    ax.set_ylabel("wave number [/mm]",fontsize=fsz);

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    bx.grid(True)
    bx.plot(bndl.freq[indx],c,"sk",markersize=8,label="phase vel.")
    bx.plot(bndl.freq[indx],cg1,"b-",markersize=8,label="group vel.")
    bx.plot(bndl.freq[indx],cg2,"r-",markersize=8,label="group vel.")

    #bx.plot(pk(k_peak[0:jmax+1]),-pk(k_peak[0:jmax+1])/k_peak[0:jmax+1],"m-")
    bx.plot(ffit,-ffit/kfit,"mo-")
    bx.set_ylim([2.5,3.5])
    bx.set_xlim([0,3])
    #bx.legend()

    Fsz=fsz+2
    ax.tick_params(labelsize=fsz)
    bx.tick_params(labelsize=fsz)
    ax.set_xlabel("frequency [MHz]",fontsize=Fsz)
    ax.set_ylabel("wave number $k_x$ [mm$^{-1}$]",fontsize=Fsz)
    bx.set_xlabel("frequency [MHz]",fontsize=Fsz)
    bx.set_ylabel("velocity [km/s]",fontsize=Fsz)

    plt.show()
    fig1.savefig("fkplot.png",bbox_inches="tight")
    fig2.savefig("vels.png",bbox_inches="tight")

