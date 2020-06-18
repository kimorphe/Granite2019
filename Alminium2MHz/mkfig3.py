import numpy as np
import matplotlib.pyplot as plt
import ascans as Asc
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar


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

    fig1=plt.figure(figsize=(14,4))    # snapshot
    fig2=plt.figure(figsize=(14,4))    # Fourier phase
    fig3=plt.figure(figsize=(14,4))    # Fourier amplitude

    fig0=plt.figure(figsize=(6,8))      # Wave Number
    fig5=plt.figure(figsize=(6,4))      # Phase Velocity

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
    """
    cnum=str(3)
    ts=[22,23,24,25]
    fs=[2.2,2.4,2.6,2.8]
    cnum=str(1)
    ts=[18,19,20,21]
    fs=[0.6,0.8,1.0,1.2]
    """
    cnum=str(4)
    ts=[26,27,28,29]
    fs=[2.2,2.4,2.6,2.8]

    cnum=str(4)
    ts=[26,27,28,29]
    fs=[3.0,3.5,4.0,5.0]
    cnum=str(2)
    ts=[20,21,22,23]
    #fs=[0.5,0.6,0.9,1.2]
    fs=[0.6,0.8,1.0,1.2]


    ax=[]; bx=[]; cx=[]; ex=[]; 
    ax=[]; bx=[]; cx=[]; ex=[]; 
    mV=1.e03
    nbin=50
    hd=["(a) ","(b) ","(c) ","(d) "];
    for k in range(4):
        ax.append(fig1.add_subplot(1,4,k+1))
        bx.append(fig2.add_subplot(1,4,k+1))
        cx.append(fig3.add_subplot(1,4,k+1))
        ax[k].tick_params(labelsize=fsz)
        bx[k].tick_params(labelsize=fsz)
        cx[k].tick_params(labelsize=fsz)
        if k>0:
            ax[k].tick_params(labelleft=False)
            bx[k].tick_params(labelleft=False)
            cx[k].tick_params(labelleft=False)

        S=cwv.get_tsnap(ts[k])  # Snapshot (vel. field)
        T=cwv.get_fsnap(fs[k])  # Spectrum (Amplitude)
        P=cwv.get_psnap(fs[k])  # Spectrum (Phase)
        Ky,Kx=grad(P)
        Kx/=-dx;
        Ky/=-dy;
        Zk=Kx+1j*Ky
        ext=[x1,x2,y2,y1]
        #   Time Snapshots
        ima=ax[k].imshow(S*mV,cmap="jet",interpolation="bilinear",extent=ext,vmin=-50,vmax=50,aspect="equal")

        #   Frequency Snapshots
        imb=bx[k].imshow(np.angle(T),aspect="equal",cmap="jet",interpolation="none",extent=ext,vmin=-np.pi,vmax=np.pi)
        txtt="t="+str(ts[k])+"[$\mu$s]"
        txtf="f="+str(fs[k])+"[MHz]"
        ax[k].set_title(hd[k]+txtt,loc="center")
        bx[k].set_title(hd[k]+txtf,loc="center")
        cx[k].set_title(hd[k]+txtf,loc="center")

        omg=fs[k]*2.*np.pi
        Amax=np.max(np.abs(T))
        imc=cx[k].imshow(np.abs(T)/Amax,aspect="equal",cmap="jet",interpolation="none",vmin=0,vmax=0.6,extent=ext)


        wgt=np.abs(Zk)
        wgt=np.reshape(wgt,[Nx*Ny])
        # K-vector PDF
        K=np.sqrt(Kx*Kx+Ky*Ky)
        hist,bins=np.histogram(np.reshape(Kx/np.pi*0.5,[Nx*Ny]),bins=nbin,range=(0.01,1.5),density=True,weights=wgt)
        binm=0.5*(bins[0:-1]+bins[1:])
        ax0.plot(binm,hist,label=txtf)
        hist,bins=np.histogram(np.reshape(np.angle(Zk)/np.pi*180,[Nx*Ny]),bins=nbin,density=True)#,weights=wgt)
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

        axdiv=make_axes_locatable(ax[k])
        bxdiv=make_axes_locatable(bx[k])
        cxdiv=make_axes_locatable(cx[k])
        cax=axdiv.append_axes("right",size="7%",pad="2%")
        cbx=bxdiv.append_axes("right",size="7%",pad="2%")
        ccx=cxdiv.append_axes("right",size="7%",pad="2%")
        cba=colorbar(ima,cax=cax)
        cbb=colorbar(imb,cax=cbx)
        cbc=colorbar(imc,cax=ccx)

        ax[k].set_xlabel("x [mm]",fontsize=fsz)
        bx[k].set_xlabel("x [mm]",fontsize=fsz)
        cx[k].set_xlabel("x [mm]",fontsize=fsz)

    #ax0.text(0,5,"(a)",fontsize=fsz+2,horizontalalignment="center")
    ax0.set_xlabel("wave number [/mm]",fontsize=fsz)
    ax0.set_ylabel("probability density",fontsize=fsz)

    #bx0.text(-170,0.012,"(b)",fontsize=fsz+2,horizontalalignment="left",verticalalignment="top")
    bx0.set_xlabel("angle [deg]",fontsize=fsz)
    bx0.set_ylabel("probability density",fontsize=fsz)
    bx0.set_xlim([-180,180])

    cx0.set_xlabel("phase velocity [km/s]",fontsize=fsz)
    cx0.set_ylabel("probability density",fontsize=fsz)

    ax0.legend()
    bx0.legend()
    cx0.legend()
    ax[0].set_ylabel("y [mm]",fontsize=12)
    bx[0].set_ylabel("y [mm]",fontsize=12)
    cx[0].set_ylabel("y [mm]",fontsize=12)

    fig0.savefig("k_pdf"+cnum+".png",bbox_inches="tight")
    fig1.savefig("snap4"+cnum+".png",bbox_inches="tight")
    fig2.savefig("phase4"+cnum+".png",bbox_inches="tight")
    fig3.savefig("ampf4"+cnum+".png",bbox_inches="tight")
    fig5.savefig("pvel"+cnum+".png",bbox_inches="tight")

    plt.show()


    #######################################################################

