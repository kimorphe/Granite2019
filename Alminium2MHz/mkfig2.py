import numpy as np
import matplotlib.pyplot as plt
import ascans as Asc


def grad(P):

    [Nx,Ny]=np.shape(P)
    PI=np.pi
    PI2=2*PI
    Px=np.zeros([Nx,Ny])
    Py=np.zeros([Nx,Ny])
    wx=np.zeros([Nx,Ny])
    wy=np.zeros([Nx,Ny])
    for i in range(Nx-1):
        for j in range(Ny):
            px=P[i+1,j]-P[i,j]
            if px > PI:
                px=-(PI2-px)
            if px < -PI:
                px=PI2+px
            Px[i,j]+=px;
            Px[i+1,j]+=px;
            wx[i,j]+=1;
            wx[i+1,j]+=1;

    for i in range(Nx):
        for j in range(Ny-1):
            py=P[i,j+1]-P[i,j]
            if py > PI:
                py=-(PI2-py)
            if py < -PI:
                py=PI2+py
            Py[i,j]+=py;
            Py[i,j+1]+=py;
            wy[i,j]+=1;
            wy[i,j+1]+=1;

    Px/=wx;
    Py/=wy;
    return(Px,Py)


if __name__=="__main__":

    fsz=12

    fig1=plt.figure(figsize=(6,6.5))
    fig2=plt.figure(figsize=(6,6.5))
    fig3=plt.figure(figsize=(6,6.5))
    fig4=plt.figure(figsize=(6,6.5))

    fig0=plt.figure(figsize=(6,8))
    fig5=plt.figure(figsize=(6,4))

    ax0=fig0.add_subplot(211)
    bx0=fig0.add_subplot(212)
    cx0=fig5.add_subplot(111)
    ax0.grid(True)
    bx0.grid(True)
    cx0.grid(True)
    ax0.tick_params(labelsize=fsz)
    bx0.tick_params(labelsize=fsz)
    cx0.tick_params(labelsize=fsz)

    dir_name="./"
    Nx=41; Ny=61
    x1=0.0; dx=0.5
    y1=15.0; dy=-0.5

    x2=x1+dx*(Nx-1)
    y2=y1+dy*(Ny-1)

    cwv=Asc.Cscan()
    cwv.load(dir_name,range(Nx*Ny),Nx,Ny)
    cwv.FFT()
    cwv.set_xaxis(x1,dx)
    cwv.set_yaxis(y1,dy)
    ts=[19,20,21,22]
    fs=[0.1,0.2,0.3,0.4]
    fs=[0.4,0.6,0.9,1.2]
    fs=[1.3,1.4,1.5,1.6]
    fs=[1.6,1.7,1.8,1.9]
    fs=[2.0,2.1,2.2,2.3]
    fs=[0.8,1.0,1.2,1.4]
    ax=[]; bx=[]; cx=[]; ex=[]; 
    mV=1.e03
    nbin=50
    hd=["(a) ","(b) ","(c) ","(d) "];
    for k in range(4):
        ax.append(fig1.add_subplot(2,2,k+1))
        bx.append(fig2.add_subplot(2,2,k+1))
        cx.append(fig3.add_subplot(2,2,k+1))
        ex.append(fig4.add_subplot(2,2,k+1))
        ax[k].tick_params(labelsize=fsz)
        bx[k].tick_params(labelsize=fsz)
        cx[k].tick_params(labelsize=fsz)
        ex[k].tick_params(labelsize=fsz)

        S=cwv.get_tsnap(ts[k])  # Snapshot (vel. field)
        T=cwv.get_fsnap(fs[k])  # Spectrum (Amplitude)
        P=cwv.get_psnap(fs[k])  # Spectrum (Phase)
        Ky,Kx=grad(P)
        Kx/=-dx;
        Ky/=-dy;
        Zk=Kx+1j*Ky
        ext=[x1,x2,y2,y1]
        #   Time Snapshots
        ima=ax[k].imshow(S*mV,cmap="jet",interpolation="bilinear",extent=ext,vmin=-10,vmax=10,aspect="equal")

        #   Frequency Snapshots
        imb=bx[k].imshow(np.angle(T),aspect="equal",cmap="jet",interpolation="none",extent=ext,vmin=-np.pi,vmax=np.pi)
        txtt="t="+str(ts[k])+"[$\mu$s]"
        txtf="f="+str(fs[k])+"[MHz]"
        ax[k].set_title(hd[k]+txtt,loc="center")
        bx[k].set_title(hd[k]+txtf,loc="center")
        cx[k].set_title(hd[k]+txtf,loc="center")
        ex[k].set_title(hd[k]+txtf,loc="center")
        """
        ax[k].text(20,14,txtt,horizontalalignment="right",verticalalignment="top",fontsize=12)
        bx[k].text(20,14,txtf,horizontalalignment="right",verticalalignment="top",fontsize=12)
        cx[k].text(20,14,txtf,horizontalalignment="right",verticalalignment="top",fontsize=12)
        ex[k].text(20,14,txtf,horizontalalignment="right",verticalalignment="top",fontsize=12)
        """

        omg=fs[k]*2.*np.pi
        imc=cx[k].imshow(omg/np.abs(Zk),aspect="equal",cmap="jet",interpolation="none",vmin=0,vmax=15.0,extent=ext)
        #ime=ex[k].imshow(np.angle(Zk)/np.pi*180.,aspect="equal",cmap="jet",interpolation="bilinear",extent=ext)
        ime=ex[k].imshow(np.angle(Zk)/np.pi*180.,aspect="equal",cmap="jet",interpolation="bilinear",extent=ext)

        wgt=np.abs(Zk)
        wgt=np.reshape(wgt,[Nx*Ny])
        # K-vector PDF
        K=np.sqrt(Kx*Kx+Ky*Ky)
        hist,bins=np.histogram(np.reshape(K/np.pi*0.5,[Nx*Ny]),bins=nbin,range=(0.01,1.5),density=True,weights=wgt)
        binm=0.5*(bins[0:-1]+bins[1:])
        ax0.plot(binm,hist,label=txtf)
        hist,bins=np.histogram(np.reshape(np.angle(Zk)/np.pi*180,[Nx*Ny]),bins=nbin,density=True,weights=wgt)
        binm=0.5*(bins[0:-1]+bins[1:])
        bx0.plot(binm,hist,label=txtf)

        # Phase Velocity PDF
        hist,bins=np.histogram(omg/np.reshape(K,[Nx*Ny]),bins=nbin,range=(0.1,8.0),density=True,weights=wgt)
        binm=0.5*(bins[0:-1]+bins[1:])
        cx0.plot(binm,hist,label=txtf)
        C=omg/K
        C=np.reshape(C,[Nx*Ny])
        Cb=np.sum(C*wgt)/np.sum(wgt[:])
        print("<c>=",np.mean(C[:]),Cb)



        fig2.colorbar(ima,ax=ax[k])
        fig2.colorbar(imb,ax=bx[k])
        fig3.colorbar(imc,ax=cx[k])
        fig4.colorbar(ime,ax=ex[k])

    ax0.text(0,5,"(a)",fontsize=fsz+2,horizontalalignment="center")
    ax0.set_xlabel("wave number [/mm]",fontsize=fsz)
    ax0.set_ylabel("probability density",fontsize=fsz)

    bx0.text(-170,0.012,"(b)",fontsize=fsz+2,horizontalalignment="left",verticalalignment="top")
    bx0.set_xlabel("angle [deg]",fontsize=fsz)
    bx0.set_ylabel("probability density",fontsize=fsz)
    bx0.set_xlim([-180,180])

    cx0.set_xlabel("phase velocity [km/s]",fontsize=fsz)
    cx0.set_ylabel("probability density",fontsize=fsz)

    ax0.legend()
    bx0.legend()
    cx0.legend()
    ax[2].set_xlabel("x [mm]",fontsize=12)
    ax[3].set_xlabel("x [mm]",fontsize=12)
    ax[0].set_ylabel("y [mm]",fontsize=12)
    ax[2].set_ylabel("y [mm]",fontsize=12)

    bx[2].set_xlabel("x [mm]",fontsize=12)
    bx[3].set_xlabel("x [mm]",fontsize=12)
    bx[0].set_ylabel("y [mm]",fontsize=12)
    bx[2].set_ylabel("y [mm]",fontsize=12)

    cx[2].set_xlabel("x [mm]",fontsize=12)
    cx[3].set_xlabel("x [mm]",fontsize=12)
    cx[0].set_ylabel("y [mm]",fontsize=12)
    cx[2].set_ylabel("y [mm]",fontsize=12)

    ex[2].set_xlabel("x [mm]",fontsize=12)
    ex[3].set_xlabel("x [mm]",fontsize=12)
    ex[0].set_ylabel("y [mm]",fontsize=12)
    ex[2].set_ylabel("y [mm]",fontsize=12)

    plt.show()

    fig0.savefig("pdf.png",bbox_inches="tight")
    fig1.savefig("snap.png",bbox_inches="tight")
    fig2.savefig("phase.png",bbox_inches="tight")

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
