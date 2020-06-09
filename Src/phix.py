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

    Kx=Img()    # create Img class instace

    num=201     # file No.
    fname="k"+str(num)+".out" # ky field data file
    fname="phix.out"
    Kx.load(fname) # load k field data


    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    y=0.5*np.arange(Kx.Ndiv[1]);
    x=15.0-0.5*np.arange(Kx.Ndiv[0]);
    ext=[y[0],y[-1],x[0],x[-1]]

    im=ax.imshow(Kx.A,extent=[-15,15,0,20],cmap="jet",interpolation="none",origin="lower")
    plt.colorbar(im)
    plt.show()
