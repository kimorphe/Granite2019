import numpy as np
import matplotlib.pyplot as plt

class Ascan:
    def blank(self,ndat):
        self.Nt=ndat
        self.t0=0.0
        self.amp=np.zeros(ndat)
        self.time=np.arange(ndat)
        self.dt=1
    def set_taxis(self,t1,t2,ndat=0):
        if ndat != 0:
            self.Nt=ndat
        self.time=np.linspace(t1,t2,self.Nt)
        self.dt=self.time[1]-self.time[0]
    def load(self,fname):
        self.fname=fname;
        fp=open(fname,"r")
        time=[]
        amp=[]
        for row in fp:
            dat=row.strip().split(",")
            time.append(float(dat[0]))
            amp.append(float(dat[1]))

        self.time=np.array(time)*1.e6
        self.amp=np.array(amp)
        self.amp-=np.mean(self.amp)

        self.dt=self.time[1]-self.time[0];
        self.Nt=len(time);

        fp.close()
    def show_t(self,ax,lwd=1.0,name=""):
        ax.plot(self.time,self.amp,linewidth=lwd,label=name)
        ax.set_xlim([self.time[0],self.time[-1]])
    def show_fft(self,ax="",lwd=1.0):
        Amp=np.fft.fft(self.amp);
        df=1/(self.dt*self.Nt)
        freq=df*np.arange(self.Nt)
        Amp[0]=0
        #ax.plot(freq,np.abs(Amp))
        Amax=np.max(np.abs(Amp))
        if ax!="":
            ax.semilogy(freq,np.abs(Amp)/Amax,linewidth=lwd)
        self.Amp=Amp
        self.freq=freq
        self.df=df
    def fft(self):
        self.Amp=np.fft.fft(self.amp);
        self.Amp[0]=0.0
        self.df=1/(self.dt*self.Nt)
        self.freq=self.df*np.arange(self.Nt)
    def show_dphi(self,ax="",f1=0.2,f2=1.5,lwd=1.0,eps=1.e-02):
        phi=np.angle(self.Amp)
        self.Weiner(eps,apply=False)
        phi*=self.Wf
        self.phi=np.unwrap(phi)
        pi2=2.*np.pi
        dw=pi2*self.df
        omg=pi2*self.freq
        nf1=np.argmin(np.abs((f1-self.freq)))
        nf2=np.argmin(np.abs((f2-self.freq)))
        tg,t0=np.polyfit(omg[nf1:nf2],self.phi[nf1:nf2],1,w=self.Wf[nf1:nf2])
        phi0=-(omg*tg+t0); 
        """
        t2,tg,t0=np.polyfit(omg[nf1:nf2],self.phi[nf1:nf2],2,w=self.Wf[nf1:nf2])
        phi0=-(omg*tg+t0+omg*omg*t2)
        """
        self.dphi=self.phi+phi0

        delt=self.dphi/(omg+1.e-05)
        #ax.plot(self.freq,self.dphi,linewidth=lwd)
        if ax!="":
            ax.plot(self.freq,delt,linewidth=lwd)
        return(tg)
    def Weiner(self,eps,apply=True,ax="",lwd=1):
        Amax=np.max(np.abs(self.Amp))
        Nw=Amax*eps;
        Pw=np.abs(self.Amp*self.Amp)
        self.Wf=Pw/(Pw+Nw*Nw)
        if apply:
            self.Amp*=self.Wf
            if ax !="":
                ax.semilogy(self.freq,np.abs(self.Amp)/Amax,linewidth=lwd)
    def show_phi(self,ax,t0=0.0,lwd=1.0,eps=1.e-02):
        phi=np.angle(self.Amp)
        self.Weiner(eps,apply=False)
        phi*=self.Wf
        self.phi=np.unwrap(phi)
        dw=2.*np.pi*self.df
        omg=2.*np.pi*self.freq
        phi0=omg*t0
        tphi=-(self.phi+phi0)/(omg+1.e-05)
        #ax.plot(self.freq,self.phi+phi0,linewidth=lwd)
        ax.plot(self.freq,tphi,linewidth=lwd)
    def get_tdly(self,eps=1.e-02,ax="",lwd=1.0):
        Amp=np.fft.fft(self.amp);
        df=1/(self.dt*self.Nt)
        freq=df*np.arange(self.Nt)
        Amp[0]=0
        self.Amp=Amp
        self.freq=freq
        self.df=df
        phi=np.angle(self.Amp)
        self.Weiner(eps,apply=False)
        phi*=self.Wf
        self.phi=np.unwrap(phi)
        #self.phi=np.unwrap(np.angle(self.Amp))
        dw=2.*np.pi*self.df
        self.tdly=-np.diff(self.phi)/dw
        self.freq2=(self.freq[0:self.Nt-1]+self.freq[1:self.Nt])*0.5
        if ax!="":
            ax.plot(self.freq2,self.tdly,linewidth=lwd)
    def LPF(self,fmax):
        Amp=np.fft.fft(self.amp);
        df=1/(self.dt*self.Nt)
        freq=df*np.arange(self.Nt)
        Amp[0]=0
        self.freq=freq;
        nfmax=int(fmax/df)
        Amp[nfmax:-1]=0.0
        Amp[0]=0;
        self.Amp=Amp
        self.df=df
        self.amp=2*np.real(np.fft.ifft(self.Amp))
    def tmax(self):
        #imax=np.argmax(abs(self.amp))
        imax=np.argmax(-self.amp)
        return(self.time[imax],self.amp[imax])
    def Butterworth(self,tb,w_6dB,mexp=4):
        tt=self.time-tb;
        Phi=1+(tt/w_6dB)**mexp
        self.amp/=Phi;
    def get_amp(self,tt):
        ii=np.argmin(np.abs(tt-self.time))
        return(self.amp[ii])

class Bscan:
    def load(self,dir_name,nums):
        awv=Ascan()
        V=np.array([])
        Nx=len(nums)
        for k in nums:
            fname=dir_name+"/scope_"+str(k)+".csv"
            awv.load(fname)
            #awv.fft()
            V=np.hstack([V,awv.amp])

        V=np.reshape(V,[Nx,awv.Nt])
        tmp=np.shape(V)
        Nx=tmp[0]; Nt=tmp[1]; 
        self.V=V
        self.Nx=Nx
        self.Nt=Nt
        self.t1=awv.time[0]
        self.t2=awv.time[-1]
        self.dt=awv.time[1]-awv.time[0]
        self.time=awv.time
        self.df=1./(self.dt*self.Nt)

        self.x1=0;
        self.x2=1;
        self.dx=0
        if Nx>1:
            self.dx=(self.x2-self.x1)/(Nx-1)
    def blank(self):
        self.V=np.array([])
        self.Nx=0
        self.Nt=0
    def add_ascan(self,awv): 
        self.Nx+=1
        self.Nt=awv.Nt
        self.t1=awv.time[0]
        self.t2=awv.time[-1]
        self.dt=awv.time[1]-awv.time[0]
        self.df=1./(self.dt*self.Nt)
        self.V=np.hstack([self.V,awv.amp])
    def finalize(self):
        self.V=np.reshape(self.V,[self.Nx,self.Nt])
        self.x1=0;
        self.x2=1;
        self.dx=0
        if self.Nx>1:
            self.dx=(self.x2-self.x1)/(self.Nx-1)

    def show(self,ax):
        Vmax=np.max(np.abs(self.V))
        t1=self.t1;
        t2=self.t2
        x1=self.x1;
        x2=self.x2;
        ext=[t1,t2,x1,x2]
        im=ax.imshow(self.V/Vmax,cmap="jet",interpolation="bilinear",extent=ext,origin="lower",vmin=-0.5,vmax=0.5,aspect="auto")
        #im=ax.imshow(self.V/Vmax,cmap="jet",interpolation="none",extent=ext,origin="lower",aspect="auto")
        self.Vmax=Vmax
        return(im)
    def show_fft(self,ax):
        print(np.shape(self.V))
        W=np.fft.fft(self.V,axis=1)
        Wmax=np.max(np.abs(W))
        W/=Wmax;
        ext=[0,self.df*(self.Nt-1),self.x1,self.x2]
        im=ax.imshow(np.log10(np.abs(W)),cmap="jet",interpolation="bicubic",extent=ext,origin="lower",aspect="auto",vmin=-4,vmax=1)
        self.W=W
        self.Wmax=Wmax
        return(im)

    def set_xaxis(self,x1,dx,npnt=0):
        self.x1=x1
        self.dx=dx
        if npnt != 0:
            self.Nx=npnt
        self.x2=self.x1+self.dx*(self.Nx-1)
        self.xcod=self.x1+np.arange(self.Nx)*self.dx
    def get_mean(self):
        awv=Ascan()
        awv.blank(self.Nt)
        awv.set_taxis(self.t1,self.t2,self.Nt)
        bsum=np.sum(self.V,0)/self.Nx
        awv.amp=bsum
        awv.fft()
        self.awv_mean=awv
        return(awv)
    def get_ascan(self,num):
        awv=Ascan()
        awv.blank(self.Nt)
        awv.set_taxis(self.t1,self.t2,self.Nt)
        awv.amp=self.V[num,:]
        return(awv)
    def normalize(self):
        Vmax=np.max(np.abs(self.V))
        self.V/=Vmax;
    def get_amp(self,tt):
        ii=np.argmin(np.abs(tt-self.time))
        return(self.V[:,ii])

class Cscan:
    def load(self,dir_name,nums,Nx,Ny):
        awv=Ascan()
        Nfile=len(nums)
        if Nfile != Nx*Ny:
            print("Inconsistent grid size !!")
            print("nfile=",Nfile)
            print("Nx x Ny= ",Nx*Ny)
            exit()

        isum=0
        for k in nums:
            fname=dir_name+"/scope_"+str(k)+".csv"
            print(fname)
            awv.load(fname)
            if isum==0:
                Nt=awv.Nt
                V=np.zeros(Nx*Ny*Nt)
            #V=np.hstack([V,awv.amp])
            i1=isum*Nt;
            i2=i1+Nt
            V[i1:i2]=awv.amp[:]
            isum+=1

        V=np.reshape(V,[Nx,Ny,Nt])
        print(np.shape(V))
        self.V=V
        self.Nx=Nx
        self.Ny=Nx
        self.Nt=Nt
        self.t1=awv.time[0]
        self.t2=awv.time[-1]
        self.dt=awv.time[1]-awv.time[0]
        self.time=awv.time
        self.df=1./(self.dt*self.Nt)
        self.freq=self.df*np.arange(self.Nt)
    def FFT(self):
        self.Z=np.fft.fft(self.V)
        self.P=np.angle(self.Z)
    def set_xaxis(self,x1,dx,npnt=0):
        self.x1=x1
        self.dx=dx
        if npnt != 0:
            self.Nx=npnt
        self.x2=self.x1+self.dx*(self.Nx-1)
        self.xcod=self.x1+np.arange(self.Nx)*self.dx
    def set_yaxis(self,y1,dy,npnt=0):
        self.y1=y1
        self.dy=dy
        if npnt != 0:
            self.Ny=npnt
        self.y2=self.y1+self.dy*(self.Ny-1)
        self.ycod=self.y1+np.arange(self.Ny)*self.dy
    def get_tsnap(self,tt):
        ii=np.argmin(np.abs(tt-self.time))
        return(np.transpose(self.V[:,:,ii]))
        #return(self.V[:,:,ii])
    def get_fsnap(self,ff):
        ii=np.argmin(np.abs(ff-self.freq))
        return(np.transpose(self.Z[:,:,ii]))
        #return(self.Z[:,:,ii])
    def get_psnap(self,ff):
        ii=np.argmin(np.abs(ff-self.freq))
        return(np.transpose(self.P[:,:,ii]))
        #return(self.P[:,:,ii])
        
        
        
if __name__=="__main__":

    
    fig=plt.figure(figsize=(8,10))
    ax=fig.add_subplot(311)
    bx=fig.add_subplot(312)
    cx=fig.add_subplot(313)
    ax.grid(True)
    bx.grid(True)
    cx.grid(True)

    bx.set_xlim([0,2.5])
    cx.set_xlim([0,2.5])
    bx.set_ylim([1.e-03,1.e01])
    #cx.set_ylim([-60,60])
    cx.set_ylim([-5,5])

    awv=Ascan()
    awv0=Ascan()
    dir_names=["X30","X50","X70","X90","X110"]
    dir_names=["./"]
    X0=30; dX=20
    Nx=len(dir_names)
    Xs=X0+np.arange(Nx)*dX
    #dir_names=["./"]
    nums=np.arange(0,61,1)
    #nums=np.arange(61,122,1)
    #nums=np.arange(122,183,1)

    tb=18.0; w_6dB=1.2;
    f1=0.1; f2=1.5
    fig0=plt.figure()
    tmp=fig0.add_subplot(111)
    tmp.set_xlim([0,2.5])

    tg0=[]
    tgx=[]
    tgx_var=[]
    eps=4.e-02
    for dir_name in dir_names:
        isum=0
        tgs=[]
        for num in nums:
            fname=dir_name+"/scope_"+str(int(num))+".csv"
            awv.load(fname)
            #awv.Butterworth(tb,w_6dB,mexp=6)
            #awv.show_t(ax)
            awv.show_fft()
            #awv.show_phi(cx,t0=tb)
            tg=awv.show_dphi(f1=f1,f2=f2,eps=eps)
            tgs.append(tg)
            if isum==0:
                awv0.blank(awv.Nt)
            awv0.amp+=awv.amp
            isum+=1
        awv0.amp/=isum
        awv0.set_taxis(awv.time[0],awv.time[-1])
        awv0.show_t(ax,lwd=2,name=dir_name)
        awv0.show_fft(ax=bx,lwd=2)
        awv0.Weiner(eps,apply=False)
        tmp.plot(awv0.freq,awv0.Wf)
        tmp.grid(True)
        tg=awv0.show_dphi(ax=cx,f1=f1,f2=f2,lwd=2,eps=eps)
        tg0.append(tg)
        print("tg=",tg)
        tgs=np.array(tgs)
        print("<tg>,std<tg>=",np.mean(tgs),np.std(tgs))
        tgx.append(np.mean(tgs))
        tgx_var.append(np.std(tgs))
    ax.legend()

    fig2=plt.figure()
    ax2=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)
    ax2.grid(True)
    bx2.grid(True)
    ax2.plot(Xs,-np.array(tg0),"-ob",markersize=10)
    ax2.plot(Xs,-np.array(tgx),"-sr",markersize=10)
    bx2.plot(Xs,np.array(tgx_var),"-sr",markersize=10)

    plt.show()

