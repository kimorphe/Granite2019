#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import plot_alph as pal 

class PSI_BNDL:
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

        fp.readline()
        amp=[]
        Pmin=[];
        Pmax=[];
        Pave=[];
        Pvar=[];
        Ppnt=[];
        for row in fp:
            dat=row.strip().split(",")
            Pmin.append(float(dat[0]))
            Pmax.append(float(dat[1]))
            Pave.append(float(dat[2]))
            Pvar.append(float(dat[3]))
            Ppnt.append(float(dat[4]))
        amp=np.array(amp)
        Nx=Nd[0]
        Ny=Nd[1]
        Pmin=np.reshape(Pmin,[Nx,Ny,Nf])
        Pmax=np.reshape(Pmax,[Nx,Ny,Nf])
        Pave=np.reshape(Pave,[Nx,Ny,Nf])
        Pvar=np.reshape(Pvar,[Nx,Ny,Nf])
        Ppnt=np.reshape(Ppnt,[Nx,Ny,Nf])

        self.Pmin=Pmin
        self.Pmax=Pmax
        self.Pave=Pave
        self.Pvar=Pvar
        self.Ppnt=Ppnt
        self.Nd=Nd;
        self.Nx=Nx
        self.Ny=Ny
        self.Nf=Nf
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

    fsz=16

    fig=plt.figure();
    ax=fig.add_subplot(111)

    V1=0; V2=200;
    fname="phix.out"
    Psi=PSI_BNDL()
    Psi.load(fname)

    freq=1.0 
    num=Psi.get_index(freq,2);
    freq=Psi.get_cod(num,2);
    print("freq=",freq,num)

    omg=freq*2.*np.pi;
    P=Psi.Pmin[:,:,num]/omg
    xx=Psi.xcod;
    yy=Psi.ycod;
    ff=Psi.freq;

    Pm=Psi.Pmin[:,:,num]/omg
    """
    Y=[]
    Pd=[]
    for i in range(Psi.Nx):
        for j in range(Psi.Ny):
            if P[i,j]<0:
                Pm[i,j]*=-10000
                continue
            Pd.append(P[i,j])
            Y.append(yy[j])
    fig0=plt.figure()
    ax0=fig0.add_subplot(111)
    ax0.plot(Y,Pd,".")
    pmin=np.min(Pm,axis=0)
    ax0.plot(yy,pmin,"og")
    ax0.grid(True)
    deg=1
    coef=np.polyfit(-yy,pmin,deg)
    coefd=np.polyder(coef)
    P_tof=np.poly1d(coef)
    Pd_tof=np.poly1d(coefd)
    #print("c=",Pd_tof(-yy));
    ax0.plot(yy,P_tof(-yy),"m")
    #print("coef=",coef);
    c=1/coef[0];
    #print("c=",c)
    c=np.sum(-yy*pmin)/np.sum(yy*yy)
    c=1/c
    #print("c=",c)
    """


    #ext=[xx[0],xx[-1],yy[0],yy[-1]]
    ext=[-yy[0],-yy[-1],-xx[0],-xx[-1]]
    im=ax.imshow(P,aspect="equal",cmap="jet",origin="lower",extent=ext,interpolation="none")

    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(im,cax=cax,orientation="vertical");
    cbar.ax.tick_params(labelsize=14)

    fname="kvec.out"
    KX=pal.BNDL()
    KX.load(fname)
    num=KX.get_index(freq,2);
    freq=KX.get_cod(num,2);
    print("freq=",freq,num)

    Kx=KX.amp[:,:,num]
    y=KX.Xa[0]+KX.dx[0]*np.arange(KX.Nd[0]);
    x=KX.Xa[1]+KX.dx[1]*np.arange(KX.Nd[1]);
    x=-x; y=-y;
    ext=[x[0],x[-1],y[0],y[-1]]
    V=np.abs(Kx);
    Kx=np.imag(Kx)+1j*np.real(Kx)
    Kx=-Kx/V;
    [X,Y]=np.meshgrid(x,y)
    X=np.transpose(X)
    Y=np.transpose(Y)
    Kx=np.transpose(Kx)
    C=np.angle(Kx)/np.pi*180.
    #ax.quiver(X+20,Y,np.real(Kx),np.imag(Kx),np.abs(C),color="k")#cmap="jet")

    ax.tick_params(labelsize=fsz)
    ax.set_xlabel("x [mm]",fontsize=fsz)
    ax.set_ylabel("y [mm]",fontsize=fsz)
    ax.set_aspect(1.0)

    """
    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    omg=Psi.freq*2.*np.pi
    [X,Omg]=np.meshgrid(np.zeros(Psi.Nx),omg)
    print(np.shape(Omg))
    print(Omg)
    V1=0
    V2=10
    jm=bx.imshow(np.transpose(Psi.Pmin[:,-2,:])/Omg,aspect="auto",extent=[xx[0],xx[-1],ff[0],ff[-1]],vmin=V1,vmax=V2,cmap="jet",origin="lower")
    bx_div=make_axes_locatable(bx);
    cbx=bx_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(jm,cax=cbx,orientation="vertical");
    """

    plt.show()
    fig.savefig("unwrapped.png",bbox_inches="tight")
