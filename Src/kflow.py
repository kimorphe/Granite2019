import numpy as np
import matplotlib.pyplot as plt
import plot_alph as pal 

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
        self.A=np.reshape(A,[Nx,Ny])
        self.xcod=Xa[0]+np.arange(Nx)*dx[0];
        self.ycod=Xa[1]+np.arange(Ny)*dx[1];
    def show(self,ax):
        ax.imshow(self.A,cmap="jet",origin="lower",aspect="equal");

class Strm:
    def load(self,fname):
        fp=open(fname,"r")

        Np=0
        ndat=[]
        nx=0
        xcod=[]
        ycod=[]
        for row in fp:
            dat=row.strip()
            if dat=="#":
                ndat.append(nx)
                nx=0
                Np+=1;
            else:
                dat=dat.split(" ")
                xcod.append(float(dat[0]))
                ycod.append(float(dat[1]))
                nx+=1

        self.Np=Np; #number of tracer particles
        self.ndat=np.array(ndat) # number of streamline data points 
        self.xcod=np.array(xcod)
        self.ycod=np.array(ycod)

    def plot(self,ax):
        
        n1=0
        for n in self.ndat:
            n2=n1+n
            xx=self.xcod[n1:n2]
            yy=self.ycod[n1:n2]
            ax.plot(xx,yy,linewidth=0.5)
            n1=n2


                

if __name__=="__main__":

    st=Strm()
    st.load("log.txt")

    Kx=Img()    # create Img class instace
    Ky=Img()    # create Img class instance

    #num=90     # file No.
    #fname1="kx"+str(num)+".out" # kx field data file
    #fname2="ky"+str(num)+".out" # ky field data file

    fname="kvec.out"
    KX=pal.BNDL()
    KX.load(fname)
    freq=0.9;
    freq=0.986366;
    num=KX.get_index(freq,2);
    freq=KX.get_cod(num,2);
    print("freq=",freq,num)
    Kx=KX.amp[:,:,num]

    #Kx.load(fname1) # load kx data
    #Ky.load(fname2) # load ky data
    C=np.abs(np.angle(-np.imag(Kx)-1j*np.real(Kx)))
    V=np.abs(Kx);
    #Kx=np.imag(Kx)+1j*np.real(Kx)
    Kx=np.real(Kx)+1j*np.imag(Kx)
    #Kx=-Kx/V;

    ## -------------- Scatter Plot ---------------
    fig3=plt.figure()
    ex=fig3.add_subplot(111)
    #[X,Y]=np.meshgrid(x,y)
    """
    K=np.abs(Kx.A+1j*Ky.A);
    Th=np.angle((Kx.A+1j*Ky.A));
    Fx=np.cos(Th); 
    Fy=np.sin(Th);
    """
    #ex.quiver(y,x,-Ky.A,-Kx.A,K,cmap="jet")
    #ex.quiver(Kx.xcod,Kx.ycod,np.transpose(Kx.A),np.transpose(Ky.A),K,cmap="jet")
    Fx=-np.real(Kx)/V;
    Fy=-np.imag(Kx)/V;
    ex.set_xlim([15,-15])
    ex.set_ylim([0,-20])
    ex.quiver(KX.xcod,KX.ycod,np.transpose(Fx),np.transpose(Fy),np.transpose(C),cmap="jet")
    qx=np.mean(Fx,axis=0)
    qy=np.mean(Fy,axis=0)
    #ex.imshow(K,extent=ext,cmap="gray",interpolation="bilinear",origin="lower",vmin=0,vmax=6)
    #ex.set_title("f="+str(Ky.freq)+"[MHz]")
    #ex.set_xlim([0,20])
    #ex.set_ylim([-15,15])
    ex.grid(True)
    ex.set_aspect(1.0)

    st.plot(ex)

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
