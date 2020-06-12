import numpy as np
import matplotlib.pyplot as plt
import ascans as Asc


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
    def get_index(self,val,axis):
        if axis==0:
            indx=np.argmin(np.abs(self.xcod-val))
        if axis==1:
            indx=np.argmin(np.abs(self.ycod-val))
        if axis==2:
            indx=np.argmin(np.abs(self.time-val))
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
            indx=indx%self.Nt;
            cod=self.time[indx]
        return(cod)
    def get_domain(self,axis):
        x=-self.xcod
        y=-self.ycod
        t= self.time
        if axis==0:
            ext=[t[0],t[-1],y[0],y[-1]]
        if axis==1:
            ext=[t[0],t[-1],y[0],y[-1]]
        if axis==2:
            ext=[y[0],y[-1],x[0],x[-1]]
        return(ext)
    def stack(self,ax):
        S=np.mean(self.amp,axis=ax)
        return(S)


if __name__=="__main__":

    dir_name="./Scopes"

    fname="scopes.csv"
    fname=dir_name+"/"+fname
    bndl=BNDL()
    bndl.load(fname)

    ext=bndl.get_domain(2);

    fsz=12
    fig1=plt.figure(figsize=(6,6.5))

    ts=[20,21,22,23]
    mV=1.e03
    hd=["(a) ","(b) ","(c) ","(d) "];
    ax=[]; 
    for k in range(4):
        ax.append(fig1.add_subplot(2,2,k+1))
        ax[k].tick_params(labelsize=fsz)
        indx=bndl.get_index(ts[k],2);
        tt=bndl.get_cod(indx,2);
        Axy=bndl.amp[:,:,indx];
        #   Time Snapshots
        ima=ax[k].imshow(Axy*mV,cmap="jet",interpolation="bilinear",extent=ext,vmin=-10,vmax=10,aspect="equal")
        fig1.colorbar(ima,ax=ax[k])

    ax[2].set_xlabel("x [mm]",fontsize=12)
    ax[3].set_xlabel("x [mm]",fontsize=12)
    ax[0].set_ylabel("y [mm]",fontsize=12)
    ax[2].set_ylabel("y [mm]",fontsize=12)


    fig2=plt.figure()
    fig3=plt.figure()
    bx=[]
    bx.append(fig2.add_subplot(111))
    bx.append(fig3.add_subplot(111))
    bx[0].tick_params(labelsize=fsz)
    bx[1].tick_params(labelsize=fsz)

    # Mean B-scan
    S=bndl.stack(0)
    ext2=bndl.get_domain(0)
    bx[0].imshow(S*mV,extent=ext2,interpolation="bilinear",aspect="auto",cmap="jet",vmin=-2,vmax=2,origin="lower")

    # B-scan at xcod 
    xcod=0
    indx=bndl.get_index(xcod,0);
    xcod=bndl.get_cod(indx,0)
    bx[1].imshow(bndl.amp[indx,:,:]*mV,extent=ext2,interpolation="bilinear",aspect="auto",cmap="jet",vmin=-4,vmax=4,origin="lower")
    bx[0].set_xlim([0,100])
    bx[1].set_xlim([0,100])
    plt.show()
    fig1.savefig("snap.png",bbox_inches="tight")
    fig2.savefig("bscan_mean.png",bbox_inches="tight")
    fig2.savefig("bscan.png",bbox_inches="tight")


    #######################################################################

    tlim=[10,90]
    nums=np.array(range(Ny))
    bwv0=Asc.Bscan()
    bwv0.blank()
    amax=[]
    amin=[]
    snap=np.array([])
    for ix in range(Nx):
        print("ix=",ix)
        bwv=Asc.Bscan()
        bwv.load(dir_name,nums)
        bwv.set_xaxis(x1,dx,npnt=Ny)
        nums+=Ny
        #-------------- A-scans -------------
        #bwv.normalize()
        awv0=bwv.get_mean()
        bwv0.add_ascan(awv0)
        amax.append(np.max(awv0.amp))
        amin.append(np.min(awv0.amp))
        snap=np.hstack([snap,bwv.get_amp(23)])
    bwv0.finalize()
    bwv0.set_xaxis(y1,dy,Ny)

    snap=np.reshape(snap,[Nx,Ny])

    bx.plot(amax,"-o",markersize=10)
    bx.plot(amin,"-s",markersize=10)


    fig3=plt.figure()
    cx=fig3.add_subplot(111)
    cx.grid(True)

    bwv0.show_fft(cx)
    im=bwv0.show(ax)
    ax.set_xlim(tlim)
    ax.tick_params(labelsize=fsz)
    ax.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    ax.set_ylabel("y [mm]",fontsize=fsz)


    fig4=plt.figure()
    ex=fig4.add_subplot(111)
    ex.imshow(snap,aspect="auto",cmap="jet",interpolation="bicubic")

    """
    fig1.savefig("bscan.png",bbox_inches="tight")
    fig2.savefig("ascan.png",bbox_inches="tight")
    """
    plt.show()
