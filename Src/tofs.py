import numpy as np
import matplotlib.pyplot as plt

class Hist:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        dat=fp.readline()
        yrng=list(map(float,dat.strip().split(",")))
        print("yrng=",yrng)

        fp.readline()
        dat=fp.readline()
        frng=list(map(float,dat.strip().split(",")))
        print("frng=",frng)

        fp.readline()
        dat=fp.readline()
        trng=list(map(float,dat.strip().split(",")))
        print("trng=",trng)

        fp.readline()
        dat=fp.readline()
        Nd=list(map(int,dat.strip().split(",")))
        print("Nd=",Nd)

        fp.readline()
        amp=[]
        for row in fp:
            amp.append(float(row))
        amp=np.array(amp)
        print(np.shape(amp))
        amp=np.reshape(amp,Nd)

        self.ycod=np.linspace(yrng[0],yrng[1],Nd[0])
        self.freq=np.linspace(frng[0],frng[1],Nd[1])
        self.time=np.linspace(trng[0],trng[1],Nd[2])
        self.Ny=Nd[0];
        self.Nf=Nd[1];
        self.Nt=Nd[0];
        self.Count=amp;

if  __name__=="__main__":
    H=Hist()
    H.load("hist_ywt.dat")

    D=np.sum(H.Count,axis=1)
    indx=np.argmax(D,axis=1)
    print(np.shape(D))
    print(np.shape(indx))


    fig=plt.figure()
    ax=fig.add_subplot(111)
    ext=[H.time[0],H.time[-1],H.ycod[0],H.ycod[-1]]
    ax.imshow(D,aspect="auto",cmap="jet",origin="lower",extent=ext,vmin=0,vmax=2000,interpolation="bilinear");
    ax.plot(H.time[indx],H.ycod,"w.")

    plt.show()
