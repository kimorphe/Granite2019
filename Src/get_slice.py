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
        fp.close()

if __name__=="__main__":
    #dir_name="../Bar2/1MHz_30x20"
    dir_name="../Alminium2MHz"
    fname="scopes.csv"

    fname=dir_name+"/"+fname

    bndl=BNDL()

    bndl.load(fname)

    #A=bndl.amp[:,10,:]
    A=bndl.amp[30,:,:]
    print(np.shape(A))

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.imshow(A,aspect="auto",cmap="jet")
    plt.show()


