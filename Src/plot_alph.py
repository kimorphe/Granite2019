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

        A=[];B=[];
        for row in fp:
            dat=row.strip().split(",")
            A.append(float(dat[0]))
            B.append(float(dat[1]))
        A=np.array(A)
        B=np.array(B)
        fp.close()

        self.Ndiv=Ndiv
        self.Xa=Xa
        self.dx=dx
        self.Wd=self.dx*self.Ndiv
        self.Xb=self.Xa+self.Wd
        self.A=np.transpose(np.reshape(A,[Nx,Ny]))
        self.B=np.transpose(np.reshape(B,[Nx,Ny]))
        self.Nx=Nx;
        self.Ny=Ny;

if __name__=="__main__":
    fname="tmp2.dat"
    alph=Img();
    alph.load(fname);
    #alph.A=np.fft.fftshift(alph.A)
    Xa=alph.Xa; 
    Ndiv=alph.Ndiv;
    dx=alph.dx;
    y=Xa[0]+dx[0]*np.arange(Ndiv[0]);
    x=Xa[1]+dx[1]*np.arange(Ndiv[1]);
    x=-x; y=-y;
    [X,Y]=np.meshgrid(x,y)
    X=np.transpose(X)
    Y=np.transpose(Y)
    Alph=alph.A/180.*np.pi


    fig=plt.figure();
    ax=fig.add_subplot(211)
    bx=fig.add_subplot(212)

    V1=-0; V2=180;

    Kx=Img();
    Kx.load("k67.out")
    Phi=np.sum(-Kx.B,axis=0);
    print(np.shape(Phi))
    y=Kx.Xa[0]+Kx.dx[0]*np.arange(Kx.Ndiv[0]);
    x=Kx.Xa[1]+Kx.dx[1]*np.arange(Kx.Ndiv[1]);
    x=-x; y=-y;
    ext=[x[0],x[-1],y[0],y[-1]]
    V=np.sqrt(Kx.A*Kx.A+Kx.B*Kx.B)
    Kx.A/=V; Kx.B/=V;
    Z=(-Kx.B-1j*Kx.A)
    #jm=bx.imshow(alph.B,aspect="equal",cmap="jet",origin="lower",extent=ext,vmin=-180.,vmax=180.,interpolation="bilinear")
    im=ax.imshow(np.abs(np.transpose(alph.A)),aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear")#,vmin=V1,vmax=V2,interpolation="none")
    ax.quiver(X+20,Y,np.cos(Alph),np.sin(Alph),np.abs(Alph/np.pi*180.),cmap="jet")
    jm=bx.imshow(np.transpose(np.abs(np.angle(Z)))/np.pi*180.,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear") #,vmin=0.,vmax=180.,interpolation="none")

    [X,Y]=np.meshgrid(x,y)
    X=np.transpose(X)
    Y=np.transpose(Y)
    A=-Kx.A; B=-Kx.B
    C=np.angle(Z)/np.pi*180.
    bx.quiver(X+20,Y,B,A,np.abs(C),cmap="jet")
    #bx.quiver(X+20,Y,B,A,cmap="jet")

    ax_div=make_axes_locatable(ax);
    bx_div=make_axes_locatable(bx);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    bax=bx_div.append_axes("right",size="5%",pad="2.5%");
    cbar1=colorbar(im,cax=cax,orientation="vertical");
    cbar2=colorbar(jm,cax=bax,orientation="vertical");
    cax.xaxis.set_ticks_position("top")
    #ax.set_xlabel("x[mm]");
    #ax.set_ylabel("y[mm]");
    #ax.set_title(fname);

    fig2=plt.figure()
    ex=fig2.add_subplot(111)
    ex.plot(y,Phi*0.5)
    ex.grid(True)
    plt.show()
    

