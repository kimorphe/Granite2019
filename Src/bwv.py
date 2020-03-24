#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

if __name__=="__main__":
    fname="bwv.out";
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


    amp=fp.read();
    amp=amp.strip();
    amp=np.array(amp.split("\n"));
    amp=amp.astype(np.float)

    ycod=Xa[0]+np.arange(Nx)*dx[0];
    tcod=Xa[1]+np.arange(Ny)*dx[1];


    ave=np.mean(amp);
    stdv=np.std(amp);
    nsig=6;
    V1=ave-nsig*stdv;
    V2=ave+nsig*stdv;
    amp=np.reshape(amp,[Nx,Ny])
    #amp=amp.transpose();



    fig,ax=plt.subplots(1,1)
    ext=[tcod[0],tcod[-1],ycod[0],ycod[-1]]
    #ax.imshow(amp,aspect="equal",cmap="jet",vmin=-0.3,vmax=0.3,origin="lower",extent=ext,interpolation="bilinear")
    im=ax.imshow(amp,aspect="auto",cmap="jet",origin="lower",extent=ext,vmin=V1,vmax=V2)
    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(im,cax=cax,orientation="vertical");
    cax.xaxis.set_ticks_position("top")
    ax.set_xlabel("time [$\mu$sec ]");
    ax.set_ylabel("y [mm]");
    ax.set_title(fname);
    plt.show()
    

