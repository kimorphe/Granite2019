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
        amp=np.reshape(amp,[Nx,Ny,Nf])
        self.amp=amp
        self.Nd=Nd;
        self.Nx=Nx
        self.Ny=Ny
        self.Nf=Nf
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


if __name__=="__main__":

    fig=plt.figure();
    #ax=fig.add_subplot(211)
    bx=fig.add_subplot(111)

    V1=-0; V2=180;

    #Kx=Img();
    #Kx.load("k100.out")
    dir_name="./"
    fname="kvec.out"
    fname=dir_name+"/"+fname
    KX=BNDL()
    KX.load(fname)
    freq=0.8;
    freq=0.9776;
    num=KX.get_index(freq,2);
    freq=KX.get_cod(num,2);
    print("freq=",freq,num)
    Kx=KX.amp[:,:,num]
    y=KX.Xa[0]+KX.dx[0]*np.arange(KX.Nd[0]);
    x=KX.Xa[1]+KX.dx[1]*np.arange(KX.Nd[1]);
    x=-x; y=-y;
    ext=[x[0],x[-1],y[0],y[-1]]
    V=np.abs(Kx);
    Kx=np.imag(Kx)+1j*np.real(Kx)
    Kx=-Kx/V;
    #im=ax.imshow(np.abs(np.transpose(alph.A)),aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="bilinear")#,vmin=V1,vmax=V2,interpolation="none")
    #ax.quiver(X+20,Y,np.cos(Alph),np.sin(Alph),np.abs(Alph/np.pi*180.),cmap="jet")
    jm=bx.imshow(np.abs(np.angle(Kx))/np.pi*180.,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="none") #,vmin=0.,vmax=180.,interpolation="none")

    [X,Y]=np.meshgrid(x,y)
    X=np.transpose(X)
    Y=np.transpose(Y)
    #A=-Kx.A; B=-Kx.B
    Kx=np.transpose(Kx)
    C=np.angle(Kx)/np.pi*180.
    #bx.quiver(X+20,Y,np.real(Kx),np.imag(Kx),np.abs(C),cmap="jet")
    bx.quiver(X+20,Y,np.real(Kx),np.imag(Kx),np.abs(C),cmap="jet")
    #bx.quiver(X+20,Y,B,A,cmap="jet")

    #ax_div=make_axes_locatable(ax);
    bx_div=make_axes_locatable(bx);
    #cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    bax=bx_div.append_axes("right",size="5%",pad="2.5%");
    #cbar1=colorbar(im,cax=cax,orientation="vertical");
    cbar2=colorbar(jm,cax=bax,orientation="vertical");
    #cax.xaxis.set_ticks_position("top")
    #ax.set_xlabel("x[mm]");
    #ax.set_ylabel("y[mm]");
    #ax.set_title(fname);

    plt.show()
    

