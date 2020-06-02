#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

class Img:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        self.freq=float(fp.readline())
        print("f=",self.freq,"[MHz]")

        fp.readline()
        dat=fp.readline()
        dat=dat.strip().split(",")
        Nx=int(dat[0])
        Ny=int(dat[1])

        fp.readline()
        dat=fp.readline()
        dat=dat.strip().split(",")
        Ndiv=np.array([Nx,Ny])
        Xa=np.zeros(2)
        print(dat)
        Xa[0]=float(dat[0]);
        Xa[1]=float(dat[1]);

        fp.readline()
        dat=fp.readline()
        dat=dat.strip().split(",")
        dx=np.zeros(2)
        dx[0]=float(dat[0]);
        dx[1]=float(dat[1]);
        fp.readline()

        A=[];B=[];C=[];
        for row in fp:
            dat=row.strip().split(",")
            A.append(float(dat[0]))
            B.append(float(dat[1]))
            C.append(float(dat[2]))
        A=np.array(A)
        B=np.array(B)
        C=np.array(C)
        fp.close()

        self.Ndiv=Ndiv
        self.Xa=Xa
        self.dx=dx
        self.Wd=self.dx*self.Ndiv
        self.Xb=self.Xa+self.Wd
        self.A=np.transpose(np.reshape(A,[Nx,Ny]))
        self.B=np.transpose(np.reshape(B,[Nx,Ny]))
        self.C=np.transpose(np.reshape(C,[Nx,Ny]))
        self.Nx=Nx;
        self.Ny=Ny;
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");

if __name__=="__main__":
    fname="h440.out";
    fp=open(fname,"r");

    H=Img();
    H.load(fname);
    fig=plt.figure()
    ax=fig.add_subplot(131)
    bx=fig.add_subplot(132)
    cx=fig.add_subplot(133)
    ext=[-10,10,0,30]

    A=H.A; B=H.B; C=H.C;
    Kmax=np.sqrt((A-B)*(A-B)*0.25+C*C)
    Kmin=-Kmax;
    Kmax+=0.5*(A+B)
    Kmin+=0.5*(A+B)
    Kmean=np.mean(Kmax[:])
    Kmean+=np.mean(Kmin[:])
    Kmean*=0.5;
    s1=np.std(Kmax[:])
    s2=np.std(Kmin[:])
    sig=(s1+s2)*0.5;
    #V1=Kmean-6*sig;
    #V2=Kmean+6*sig;
    V1=-20;V2=20
    im0=ax.imshow(Kmin,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear",vmin=V1,vmax=V2)
    im1=bx.imshow(Kmax,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear",vmin=V1,vmax=V2)
    sgn=Kmax*Kmin/abs(Kmax*Kmin)
    im2=cx.imshow(np.sqrt(abs(Kmax*Kmin))*sgn,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear",vmin=V1,vmax=V2)
    #im2=cx.imshow(Kmax*Kmin/abs(Kmax*Kmin),aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear")

    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    #cax.xaxis.set_ticks_position("top")
    cbar=colorbar(im0,cax=cax,orientation="vertical");

    bx_div=make_axes_locatable(bx);
    cax=bx_div.append_axes("right",size="5%",pad="2.5%");
    #cax.xaxis.set_ticks_position("top")
    cbar=colorbar(im1,cax=cax,orientation="vertical");

    cx_div=make_axes_locatable(cx);
    cax=cx_div.append_axes("right",size="5%",pad="2.5%");
    #cax.xaxis.set_ticks_position("top")
    cbar=colorbar(im2,cax=cax,orientation="vertical");

    fig2=plt.figure()
    ax2=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)
    Kmax=np.reshape(Kmax,[H.Nx*H.Ny])
    Kmin=np.reshape(Kmin,[H.Nx*H.Ny])
    #ax2.hist(Kmax,bins=40)
    ax2.hist(Kmax-Kmin,bins=50,range=[-40,40])
    bx2.hist((Kmax+Kmin)*0.5,bins=50,range=[-40,40])
    ax2.grid(True)
    bx2.grid(True)

    """
    ext=[-10,10,0,30]
    #ax.imshow(amp,aspect="equal",cmap="jet",vmin=-0.3,vmax=0.3,origin="lower",extent=ext,interpolation="bilinear")
    im=ax.imshow(amp,aspect="equal",cmap="jet",origin="lower",extent=ext,vmin=V1,vmax=V2)

    ax.set_xlabel("x[mm]");
    ax.set_ylabel("y[mm]");
    ax.set_title(fname);
    """
    plt.show()
    

