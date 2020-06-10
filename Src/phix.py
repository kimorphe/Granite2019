import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import plot_alph as pal 

class PSI_BNDL:
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

        fp.readline()
        amp=[]
        Pmin=[];
        Pmax=[];
        Pave=[];
        Pvar=[];
        for row in fp:
            dat=row.strip().split(",")
            Pmin.append(float(dat[0]))
            Pmax.append(float(dat[1]))
            Pave.append(float(dat[2]))
            Pvar.append(float(dat[3]))
        amp=np.array(amp)
        Nx=Nd[0]
        Ny=Nd[1]
        Pmin=np.reshape(Pmin,[Nx,Ny,Nf])
        Pmax=np.reshape(Pmax,[Nx,Ny,Nf])
        Pave=np.reshape(Pave,[Nx,Ny,Nf])
        Pvar=np.reshape(Pvar,[Nx,Ny,Nf])

        self.Pmin=Pmin
        self.Pmax=Pmax
        self.Pave=Pave
        self.Pvar=Pvar
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
    ax=fig.add_subplot(111)

    V1=0; V2=200;
    fname="phix.out"
    Psi=PSI_BNDL()
    Psi.load(fname)

    freq=0.9776;
    freq=2.4
    num=Psi.get_index(freq,2);
    freq=Psi.get_cod(num,2);
    print("freq=",freq,num)

    omg=freq*2.*np.pi;
    #P=np.transpose(Psi.Pave[:,:,num])/omg
    P=np.transpose(Psi.Pmin[:,:,num])/omg
    xx=Psi.xcod;
    yy=Psi.ycod;
    ff=Psi.freq;
    ext=[xx[0],xx[-1],yy[0],yy[-1]]
    im=ax.imshow(P,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="none")

    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(im,cax=cax,orientation="vertical");

    fname="kvec.out"
    KX=pal.BNDL()
    KX.load(fname)
    num=KX.get_index(freq,2);
    freq=KX.get_cod(num,2);
    print("freq=",freq,num)
    Kx=KX.amp[:,:,num]
    C=np.abs(np.angle(-np.imag(Kx)-1j*np.real(Kx)))
    V=np.abs(Kx);
    Fx=-np.real(Kx)/V;
    Fy=-np.imag(Kx)/V;
    #ax.quiver(KX.xcod,KX.ycod,np.transpose(Fx),np.transpose(Fy),np.transpose(C),cmap="jet")
    ax.quiver(KX.xcod,KX.ycod,np.transpose(Fx),np.transpose(Fy),color="k")

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    omg=Psi.freq*2.*np.pi
    [X,Omg]=np.meshgrid(np.zeros(Psi.Nx),omg)
    print(np.shape(Omg))
    print(Omg)
    V1=0
    V2=10
    jm=bx.imshow(np.transpose(Psi.Pmin[:,-2,:])/Omg,aspect="auto",extent=[xx[0],xx[-1],ff[0],ff[-1]],vmin=V1,vmax=V2,cmap="jet",origin="lower")
    bx_div=make_axes_locatable(bx);
    cbx=bx_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(jm,cax=cbx,orientation="vertical");

    plt.show()

"""
class ImgS:
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
            dat=row.strip(); #.split(",")
            A.append(float(dat))
            #B.append(float(dat[1]))
        A=np.array(A)
        B=np.array(B)
        fp.close()

        self.Ndiv=Ndiv
        self.Xa=Xa
        self.dx=dx
        self.Wd=self.dx*self.Ndiv
        self.Xb=self.Xa+self.Wd
        self.A=np.transpose(np.reshape(A,[Nx,Ny]))
        #self.B=np.transpose(np.reshape(B,[Nx,Ny]))
        self.Nx=Nx;
        self.Ny=Ny;
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");

if __name__=="__main__":

    P=ImgS()    # create Img class instace

    num=201     # file No.
    fname="k"+str(num)+".out" # ky field data file
    fname="phix.out"
    P.load(fname) # load k field data


    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    y=-0.5*np.arange(P.Ndiv[1]);
    x=15.0-0.5*np.arange(P.Ndiv[0]);
    ext=[x[0],x[-1],y[0],y[-1]]

    im=ax.imshow(P.A,extent=ext,cmap="jet",interpolation="none",origin="lower")
    ax.grid(True)
    plt.colorbar(im)


    qx=np.mean(Fx,axis=0)
    qy=np.mean(Fy,axis=0)
    #ex.grid(True)
    #ex.set_aspect(1.0)

    fig2=plt.figure()
    ax=fig2.add_subplot(111)
    ax.plot(KX.ycod,qx)
    ax.plot(KX.ycod,qy)
    ax.plot(KX.ycod,np.abs(qx+1j*qy))
    ax.grid(True)
    Qx=np.mean(qx)
    Qy=np.mean(qy)
    print(Qx,Qy,np.abs(Qx+1j*Qy))
    plt.show()

"""
