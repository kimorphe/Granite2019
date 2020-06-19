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
        self.A=np.transpose(np.reshape(A,[Nx,Ny]))
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");
class Stats:
    def load(self,fname):
        fp=open(fname,"r")
        print(fp.readline())
        freq=[]
        kb=[]
        sk=[];
        thb=[]
        sth=[];
        for row in fp:
            dat=row.strip().split(" ")
            freq.append(float(dat[0]))
            kb.append(float(dat[1]))
            sk.append(float(dat[2]))
            thb.append(float(dat[3]))
            sth.append(float(dat[4]))


        self.freq=np.array(freq)
        self.kb=np.array(kb)
        self.thb=np.array(thb)
        self.sk=np.array(sk)
        self.sth=np.array(sth)

        print(self.thb)
    def kwfit(self,f1,f2,deg):
        freq=self.freq;
        indx1=np.argmin(np.abs(freq-f1))
        indx2=np.argmin(np.abs(freq-f2))

        freq=self.freq[indx1:indx2+1]
        kb=self.kb[indx1:indx2+1]

        coef=np.polyfit(freq,kb,deg)
        pkw=np.poly1d(coef)
        coefd=np.polyder(coef)
        pkwd=np.poly1d(coefd)

        kfit=pkw(freq)
        kdfit=pkwd(freq)

        self.fbnd=freq
        self.kfit=kfit
        self.c=freq/kfit
        self.cg=1/kdfit





if __name__=="__main__":

    stat=Stats()
    stat.load("mean.out")

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

    kv=np.arange(K.Ndiv[1])*K.dx[1]+K.Xa[1];
    alph=np.arange(A.Ndiv[1])*A.dx[1]+A.Xa[1];
    freq=np.arange(K.Ndiv[0])*K.dx[0]+K.Xa[0];

    #[X,Y]=np.meshgrid(x,y)
    ext1=[freq[0],freq[-1],kv[0],kv[-1]]
    ext2=[freq[0],freq[-1],alph[0]+90,alph[-1]+90]
    imb=bx.imshow(A.A,extent=ext2,cmap="jet",interpolation="bilinear",origin="lower",aspect="auto")
    bx.grid(True)
    #ex.set_aspect(1.0)

    ima=ax.imshow(K.A,extent=ext1,cmap="jet",interpolation="bilinear",origin="lower",aspect="auto",vmin=0,vmax=0.3)
    ax.grid(True)
    fsz=16
    ax.set_xlabel("frequency [MHz]",fontsize=fsz)
    ax.set_ylabel("wave number [mm$^{-1}$]",fontsize=fsz)
    bx.set_xlabel("frequency [MHz]",fontsize=fsz)
    bx.set_ylabel("wave direction [deg]",fontsize=fsz)
    ax.tick_params(labelsize=fsz-2)
    bx.tick_params(labelsize=fsz-2)

    fig1.colorbar(ima)
    fig2.colorbar(imb)

    lwd=1.0
    stat.thb+=90;
    ax.plot(stat.freq,stat.kb,"w",linewidth=lwd)
    bx.plot(stat.freq,stat.thb,"w",linewidth=lwd)

    lwd=1.0
    ax.plot(stat.freq,stat.kb+stat.sk,"m",linewidth=lwd)
    ax.plot(stat.freq,stat.kb-stat.sk,"m",linewidth=lwd)
    bx.plot(stat.freq,stat.thb+stat.sth,"m-",linewidth=lwd)
    bx.plot(stat.freq,stat.thb-stat.sth,"m-",linewidth=lwd)

    #ax.plot(stat.fbnd,stat.kfit,"k")

    f1=stat.freq[0];
    f2=stat.freq[-1];
    f1=0.0; f2=3.0
    k1=0.0; k2=1.2
    ax.set_xlim([f1,f2])
    ax.set_ylim([k1,k2])
    bx.set_xlim([f1,f2])

    # Polynomial Approximation
    """
    deg=2
    stat.kwfit(0.6,1.6,deg)
    fig3=plt.figure()
    cx=fig3.add_subplot(111)
    cx.grid(True)
    cx.set_ylim([2,4])
    cx.plot(stat.fbnd,stat.c,"-",markersize=8,label="phase vel.")
    cx.plot(stat.fbnd,stat.cg,"-",markersize=8,label="group vel.")
    cx.legend()
    """

    plt.show()
    fig1.savefig("kw.png",bbox_inches="tight")
    fig2.savefig("thw.png",bbox_inches="tight")

