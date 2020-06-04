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
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");

if __name__=="__main__":

    Kx=Img()    # create Img class instace

    num=201     # file No.
    fname="k"+str(num)+".out" # ky field data file
    Kx.load(fname) # load k field data


    ## -------------- Scatter Plot ---------------
    fig1=plt.figure()
    ax=fig1.add_subplot(211)
    ax2=fig1.add_subplot(212)
    ndat=Kx.Ndiv[0]*Kx.Ndiv[1]; # Number of data
    xi0=np.reshape(-Kx.A,[ndat,1])  # straighten kx data
    xi1=np.reshape(-Kx.B,[ndat,1])  # straighten ky data
    ax.plot(xi0,xi1,".",markersize=2)   #  show scatter plot of (kx,ky)
    ax.grid(True)
    ax.set_aspect(1.0)
    ax.set_xlabel("kx[/mm]");
    ax.set_ylabel("ky[/mm]");


    fig2=plt.figure()
    bx=fig2.add_subplot(211)
    cx=fig2.add_subplot(212)
    xi=np.sqrt(xi0*xi0+xi1*xi1) # wave number (magnitude)
    alph=np.angle(xi0+1j*xi1)/np.pi*180.    # wave number direction

    ax2.plot(alph,xi,".",markersize=2)
    ax2.grid(True)
    #ax2.set_aspect(1.0)
    ax2.set_xlabel("angle")
    ax2.set_ylabel("wave number")

    nbin=50
    hist_xi,bins=np.histogram(xi,bins=nbin,range=(0.01,1.5),density=False) # weights=wgt
    bin_xi=0.5*(bins[0:-1]+bins[1:])
    #bx.hist(xi,bins=50)
    bx.plot(bin_xi,hist_xi,label=fname)
    bx.grid(True);

    hist_alph,bins=np.histogram(alph,bins=nbin,density=False) #weight=wgt
    bin_alph=(bins[0:-1]+bins[1:])*0.5
    cx.hist(alph,bins=50)
    cx.plot(bin_alph,hist_alph)
    cx.grid(True);
    bx.set_xlabel("wave number [/mm]")
    cx.set_xlabel("wave number angle [deg]")
    bx.set_ylabel("count")
    cx.set_ylabel("count")


    fig3=plt.figure()
    ex=fig3.add_subplot(111)
    y=0.5*np.arange(Kx.Ndiv[1]);
    x=15.0-0.5*np.arange(Kx.Ndiv[0]);
    print("x=",x)
    print("y=",y)
    #[X,Y]=np.meshgrid(x,y)
    C=np.sqrt(Kx.A*Kx.A+Kx.B*Kx.B)
    K=np.abs(Kx.A+1j*Kx.B);
    A=np.angle(-(Kx.A+1j*Kx.B));
    Fx=np.cos(A)/K;
    Fy=np.sin(A)/K;
    ext=[y[0],y[-1],x[0],x[-1]]

    ex.imshow(A,extent=ext,cmap="jet",interpolation="bilinear",origin="lower")
    plt.show()
    ex.quiver(y,x,-Ky.A,-Kx.A,C,cmap="jet")
    #ex.imshow(C, extent=ext, cmap="gray",origin="lower")
    ex.imshow(K,extent=ext,cmap="gray",interpolation="bilinear",origin="lower",vmin=0,vmax=6)
    #ex.quiver(y,x,Fx,Fy,1./K,cmap="jet")
    ex.set_title("f="+str(Ky.freq)+"[MHz]")
    ex.set_xlim([0,20])
    ex.set_ylim([-15,15])
    ex.grid(True)
    ex.set_aspect(1.0)

    fig4=plt.figure()
    fx=fig4.add_subplot(111)
    fx.hist(1./xi,bins=50,range=[0,10])
    fx.grid(True)
    fx.set_xlabel("phase velocity [km/s]")
    fx.set_ylabel("cell count");
    plt.show()
