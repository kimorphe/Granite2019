import numpy as np
import matplotlib.pyplot as plt


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
        t1=dat[0]
        dt=dat[1]
        Nt=int(dat[2])
        print("t1,dt,Nt=",t1,dt,Nt)

        fp.readline()
        amp=[]
        for row in fp:
            amp.append(float(row))
        amp=np.array(amp)
        Nx=Nd[0]
        Ny=Nd[1]
        #amp=amp.astype(np.float)
        amp=np.reshape(amp,[Nx,Ny,Nt])
        print(np.shape(amp))
        self.amp=amp
        self.Nx=Nx
        self.Ny=Ny
        self.Nt=Nt
        self.time=np.arange(Nt)*dt+t1;
        self.dt=dt;
        self.Xa=Xa
        self.dx=dx
        self.xcod=np.arange(Nx)*dx[0]+Xa[0]
        self.ycod=np.arange(Ny)*dx[1]+Xa[1]
        fp.close()

if __name__=="__main__":
    dir_name="../Alminium2MHz"
    fname="scopes_win.csv"

    dir_name="../Bar2/1MHz_30x20"
    dir_name="../CoreM_short3/x30y20"
    fname="scopes.csv"

    fname=dir_name+"/"+fname

    bndl=BNDL()

    bndl.load(fname)

    #A=bndl.amp[:,10,:]
    #A=bndl.amp[60,:,:]
    A=bndl.amp[:,:,460]
    print(np.shape(A))

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ext=[bndl.time[0],bndl.time[-1],bndl.ycod[-1],bndl.ycod[0]];
    ax.imshow(A,aspect="auto",cmap="jet",extent=ext)
    plt.show()


