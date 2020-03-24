#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

if __name__=="__main__":
    #fname="Linf.out";
    fname="L2.out";
    fp=open(fname,"r");
    dat=fp.readline()
    dat=list(map(int,dat.strip().split(",")))

    Nx=dat[0];
    Ny=dat[1];

    amp=fp.read();
    amp=amp.strip();
    amp=np.array(amp.split("\n"));
    amp=amp.astype(np.float)


    ave=np.mean(amp);
    stdv=np.std(amp);
    nsig=2;
    V1=ave-nsig*stdv;
    V2=ave+nsig*stdv;
    amp=np.reshape(amp,[Nx,Ny])
    amp=amp.transpose();



    fig,ax=plt.subplots(1,1)
    ext=[-10,10,0,30]
    #ax.imshow(amp,aspect="equal",cmap="jet",vmin=-0.3,vmax=0.3,origin="lower",extent=ext,interpolation="bilinear")
    im=ax.imshow(amp,aspect="equal",cmap="jet",origin="lower",extent=ext,vmin=V1,vmax=V2)
    ax_div=make_axes_locatable(ax);
    cax=ax_div.append_axes("right",size="5%",pad="2.5%");
    cbar=colorbar(im,cax=cax,orientation="vertical");
    cax.xaxis.set_ticks_position("top")
    ax.set_xlabel("x[mm]");
    ax.set_ylabel("y[mm]");
    ax.set_title(fname);
    plt.show()
    

