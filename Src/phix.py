import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import plot_alph as pal 


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

    fname="kvec.out"
    KX=pal.BNDL()
    KX.load(fname)
    freq=P.freq;
    num=KX.get_index(freq,2);
    freq=KX.get_cod(num,2);
    print("freq=",freq,num)
    Kx=KX.amp[:,:,num]
    C=np.abs(np.angle(-np.imag(Kx)-1j*np.real(Kx)))
    V=np.abs(Kx);
    Kx=np.real(Kx)+1j*np.imag(Kx)

    #fig3=plt.figure()
    #ex=fig3.add_subplot(111)
    Fx=-np.real(Kx)/V;
    Fy=-np.imag(Kx)/V;
    #ex.set_xlim([15,-15])
    #ex.set_ylim([0,-20])
    ax.quiver(KX.xcod,KX.ycod,np.transpose(Fx),np.transpose(Fy),np.transpose(C),cmap="jet")
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
