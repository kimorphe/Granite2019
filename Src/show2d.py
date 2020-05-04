import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

if __name__=="__main__":
    fname="tg.out";
    fname="Linf.out";

    fp=open(fname,"r");
    fp.readline()
    dat=fp.readline()
    dat=list(map(int,dat.strip().split(",")))
    Nx=dat[0]; Ny=dat[1];

    fp.readline()
    dat=fp.readline()
    Xa=list(map(float,dat.strip().split(",")));

    fp.readline()
    dat=fp.readline()
    dx=list(map(float,dat.strip().split(",")));
    fp.readline()

    xcod=Xa[0]+np.arange(Nx)*dx[0];
    ycod=Xa[1]+np.arange(Ny)*dx[1];

    amp=fp.read();
    amp=amp.strip();
    amp=np.array(amp.split("\n"));
    amp=amp.astype(np.float)
    amp=np.reshape(amp,[Nx,Ny])
    amp=amp.transpose();

    fig,ax=plt.subplots(1,1)
    ext=[xcod[0],xcod[-1],ycod[0],ycod[-1]]

    ave=np.mean(amp);
    stdv=np.std(amp);
    nsig=3;
    V1=ave-nsig*stdv;
    V2=ave+nsig*stdv;
    #ax.imshow(amp,aspect="equal",cmap="jet",vmin=-0.3,vmax=0.3,origin="lower",extent=ext,interpolation="bilinear")
    im=ax.imshow(amp,aspect="equal",cmap="jet",origin="lower",extent=ext)
    #im=ax.imshow(amp,aspect="equal",cmap="jet",vmin=V1,vmax=V2,origin="lower",extent=ext)
    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(im,cax=cax,orientation="vertical");
    cax.xaxis.set_ticks_position("top")
    ax.set_xlabel("x [mm]");
    ax.set_ylabel("y [mm]");
    ax.set_title(fname);

    del_tg=np.gradient(amp);

    fig2,ax2=plt.subplots(1,3);
    nsig=2;
    Tg=np.sqrt(del_tg[0]**2+del_tg[1]**2)
    for k in range(2):
        ave=np.mean(del_tg[k]);
        stdv=np.std(del_tg[k]);
        V1=ave-nsig*stdv;
        V2=ave+nsig*stdv;
        ax2[k].imshow(del_tg[k],origin="lower",cmap="jet",aspect="equal",vmin=V1,vmax=V2,extent=ext,interpolation="bicubic");

    ave=np.mean(Tg);
    stdv=np.std(Tg);
    V1=ave-nsig*stdv;
    V2=ave+nsig*stdv;
    ax2[2].imshow(Tg,origin="lower",cmap="jet",aspect="equal",extent=ext,interpolation="bilinear",vmin=V1,vmax=V2);


    fig3,bx=plt.subplots(1,1)
    nums=[0,5,10,15,20,25,30,25,40]
    nums=range(41);
    tgbx=np.mean(amp,1);
    for k in nums:
        amp[:,k]-=tgbx;

    ave=np.mean(amp);
    stdv=np.std(amp);
    nsig=2;
    V1=ave-nsig*stdv;
    V2=ave+nsig*stdv;
    bx.imshow(-amp,cmap="jet",vmin=V1,vmax=V2,origin="lower",interpolation="none",extent=ext);
    #bx.plot(tgbx,"k",linewidth=3);
    bx.grid(True)
    plt.show()
    

