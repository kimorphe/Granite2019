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
        return(cod)
    def get_domain(self,axis):
        x=-self.xcod
        y=-self.ycod
        f= self.freq
        if axis==0:
            ext=[f[0],f[-1],y[0],y[-1]]
        if axis==1:
            ext=[f[0],f[-1],y[0],y[-1]]
        if axis==2:
            ext=[y[0],y[-1],x[0],x[-1]]
        return(ext)
    def stack(self,ax):
        S=np.mean(self.amp,axis=ax)
        return(S)



if __name__=="__main__":

    dir_name="./Scopes"

    fname="scopes.fft"
    fname=dir_name+"/"+fname
    bndl=BNDL()
    bndl.load(fname)

    ext=bndl.get_domain(2);

    fsz=12
    fig1=plt.figure(figsize=(6,6.5))

    fs=[0.5, 0.7, 0.9, 1.1]
    mV=1.e03
    hd=["(a) ","(b) ","(c) ","(d) "];
    ax=[];
    #   Snapshots
    for k in range(4):
        ax.append(fig1.add_subplot(2,2,k+1))
        ax[k].tick_params(labelsize=fsz)
        axdiv=make_axes_locatable(ax[k])
        cax=axdiv.append_axes("right",size="7%",pad="2%")
        indx=bndl.get_index(fs[k],2);
        ff=bndl.get_cod(indx,2);
        Phi=np.angle(bndl.amp[:,:,indx]);
        ima=ax[k].imshow(Phi,cmap="jet",interpolation="none",extent=ext,vmin=-np.pi,vmax=np.pi,aspect="equal")
        cba=colorbar(ima,cax=cax)

        txtf="f="+str(fs[k])+"[MHz]"
        ax[k].set_title(hd[k]+txtf,loc="center")

    ax[2].set_xlabel("x [mm]",fontsize=12)
    ax[3].set_xlabel("x [mm]",fontsize=12)
    ax[0].set_ylabel("y [mm]",fontsize=12)
    ax[2].set_ylabel("y [mm]",fontsize=12)


    fig2=plt.figure(figsize=(6,3.5))
    fig3=plt.figure(figsize=(6,3.5))
    bx=[]
    bx.append(fig2.add_subplot(111))
    bx.append(fig3.add_subplot(111))


    # Mean B-scan
    S=bndl.stack(0)
    ext2=bndl.get_domain(0)
    imb0=bx[0].imshow(np.abs(S*mV),extent=ext2,interpolation="bilinear",aspect="auto",cmap="jet",origin="lower")
    bxdiv=make_axes_locatable(bx[0])
    cbx=bxdiv.append_axes("right",size="7%",pad="2%")
    cb0=colorbar(imb0,cax=cbx)

    # B-scan at xcod 
    xcod=0
    indx=bndl.get_index(xcod,0);
    xcod=bndl.get_cod(indx,0)
    imb1=bx[1].imshow(np.abs(bndl.amp[indx,:,:]*mV),extent=ext2,interpolation="bilinear",aspect="auto",cmap="jet",origin="lower")
    bxdiv=make_axes_locatable(bx[1])
    cbx=bxdiv.append_axes("right",size="7%",pad="2%")
    cb0=colorbar(imb1,cax=cbx)

    for k in range(2):
        bx[k].tick_params(labelsize=fsz)
        bx[k].set_xlim([0,2.5])
        bx[k].set_xlabel("freq [MHz]",fontsize=fsz);
        bx[k].set_ylabel("x [mm]",fontsize=fsz);
    plt.show()
    fig1.savefig("fsnap.png",bbox_inches="tight")
    fig2.savefig("spctrg_mean.png",bbox_inches="tight")
    fig2.savefig("spctrg.png",bbox_inches="tight")

