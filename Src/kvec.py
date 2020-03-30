import numpy as np
import matplotlib.pyplot as plt

class Img:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        self.freq=float(fp.readline())

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

        A=[];
        for row in fp:
            A.append(float(row))
        A=np.array(A)
        fp.close()

        self.Ndiv=Ndiv
        self.Xa=Xa
        self.dx=dx
        self.Wd=self.dx*self.Ndiv
        self.Xb=self.Xa+self.Wd
        self.A=np.reshape(A,[Nx,Ny])
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");

if __name__=="__main__":

    Kx=Img()
    Ky=Img()

    fname1="kx13.out"
    fname2="ky13.out"

    Kx.load(fname1)
    Ky.load(fname2)
    print(Kx.Wd)


    fig=plt.figure()
    ax=fig.add_subplot(111)
    #Kx.show(ax) 

    ndat=Kx.Ndiv[0]*Kx.Ndiv[1];
    xi0=np.reshape(Kx.A,[ndat,1])
    xi1=np.reshape(Ky.A,[ndat,1])
    ax.plot(xi0,xi1,".",markersize=2)
    ax.grid(True)
    ax.set_aspect(1.0)

    fig2=plt.figure()
    bx=fig2.add_subplot(211)
    cx=fig2.add_subplot(212)
    xi=np.sqrt(xi0*xi0+xi1*xi1)
    alph=np.angle(xi0+1j*xi1)
    bx.hist(xi,bins=40)
    bx.grid(True);
    cx.hist(alph,bins=40)
    cx.grid(True);


    fig3=plt.figure()
    ex=fig3.add_subplot(111)
    y=0.5*np.arange(Kx.Ndiv[1]);
    x=15.0-0.5*np.arange(Kx.Ndiv[0]);
    print("x=",x)
    print("y=",y)
    #[X,Y]=np.meshgrid(x,y)
    C=np.sqrt(Ky.A*Ky.A+Kx.A*Kx.A)
    K=np.abs(Ky.A+1j*Kx.A);
    A=np.angle(-(Ky.A+1j*Kx.A));
    Fx=np.cos(A)/K;
    Fy=np.sin(A)/K;
    #ex.quiver(y,x,-Ky.A,-Kx.A,C,cmap="jet")
    ex.quiver(y,x,Fx,Fy,K,cmap="jet")
    ex.set_title("f="+str(Ky.freq)+"[MHz]")
    ex.set_xlim([0,20])
    ex.set_ylim([-15,15])
    ex.grid(True)
    ex.set_aspect(1.0)
    plt.show()
