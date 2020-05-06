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

    K=Img()    # create Img class instace
    A=Img()    # create Img class instance

    fname1="Klen.out";
    fname2="Kalp.out";

    K.load(fname1) # load kx data
    A.load(fname2) # load ky data

    fig1=plt.figure()
    ax=fig1.add_subplot(111)

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    y=0.5*np.arange(K.Ndiv[1]);
    x=15.0-0.5*np.arange(K.Ndiv[0]);

    kv=np.arange(K.Ndiv[0])*K.dx[0]+K.Xa[0];
    alph=np.arange(A.Ndiv[0])*A.dx[0]+A.Xa[0];
    freq=np.arange(K.Ndiv[1])*K.dx[1]+K.Xa[1];

    print("x=",x)
    print("y=",y)
    #[X,Y]=np.meshgrid(x,y)
    ext1=[freq[0],freq[-1],kv[0],kv[-1]]
    ext2=[freq[0],freq[-1],alph[0],alph[-1]]
    bx.imshow(A.A,extent=ext2,cmap="jet",interpolation="bilinear",origin="lower",aspect="auto")
    bx.set_title("f="+str(K.freq)+"[MHz]")
    bx.grid(True)
    #ex.set_aspect(1.0)

    ax.imshow(K.A,extent=ext1,cmap="jet",interpolation="bilinear",origin="lower",aspect="auto")
    ax.set_title("f="+str(K.freq)+"[MHz]")
    ax.grid(True)

    plt.show()
