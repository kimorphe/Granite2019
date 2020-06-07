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
        f1=dat[0]
        df=dat[1]
        Nf=int(dat[2])
        print("f1,df,Nf=",f1,df,Nf)
        #Nf=int(Nf/2)

        fp.readline()
        amp=[]
        for row in fp:
            dat=row.strip().split(",")
            z=float(dat[0])+1j*float(dat[1])
            amp.append(z)
        amp=np.array(amp)
        Nx=Nd[0]
        Ny=Nd[1]
        #amp=amp.astype(np.float)
        print(np.size(amp))
        print(Nx*Ny*Nf)
        amp=np.reshape(amp,[Nx,Ny,Nf])
        print(np.shape(amp))
        self.amp=amp
        self.Nx=Nx
        self.Ny=Ny
        self.Nf=Nf
        #self.time=np.arange(Nt)*dt+t1;
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
    dir_name="../Alminium2MHz"
    fname="scopes_win.csv"

    dir_name="../CoreM_short3/x30y20"
    fname="scopes.fft"

    fname=dir_name+"/"+fname

    bndl=BNDL()

    print(fname)
    bndl.load(fname)
    freq=0.8;
    num=bndl.get_index(freq,2);
    freq=bndl.get_cod(num,2);

    fig=plt.figure()
    ax=fig.add_subplot(211)
    bx=fig.add_subplot(212)
    f1=bndl.freq[0];
    f2=bndl.freq[-1];
    x1=bndl.xcod[0];
    x2=bndl.xcod[-1];
    y1=bndl.ycod[0];
    y2=bndl.ycod[-1];
    ext=[x1,x2,y1,y2]
    nums=[num]
    isum=0
    [Y,X]=np.meshgrid(bndl.ycod,bndl.xcod)
    #X=np.transpose(X)
    #Y=np.transpose(Y)
    for k in nums:
        #yy=bndl.ycod[k]
        A=bndl.amp[:,:,k]
        W=np.abs(A)
        A=np.angle(A)
        if isum==0:
            Wmax=np.max(W[:])
        #ax.imshow(W,aspect="equal",cmap="jet",\
        ax.imshow(np.transpose(W),aspect="equal",cmap="jet",\
        extent=ext,origin="lower",interpolation="none",vmin=0,vmax=Wmax*0.3)
        bx.imshow(np.transpose(A),aspect="equal",cmap="jet",\
        extent=ext,origin="lower",interpolation="none")
        #ax.quiver(X,Y,-np.real(A),-np.imag(A))#,np.angle(A),cmap="jet")
        #txt="y="+str(yy)+"[mm]" ax.set_title(txt)
        fnout="f"+str(isum)+"png"
        #ax.set_ylim([0,1.6])
        #fig.savefig(fnout,bbox_inches="tight")
        #plt.cla()
        print(np.shape(X))
        print(np.shape(A))
        ax.grid(True)
        bx.grid(True)
        isum+=1;
    plt.show()


