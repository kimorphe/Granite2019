import numpy as np
import matplotlib.pyplot as plt
import ascans as Asc


if __name__=="__main__":


    #-------------- B-scan image -------------
    dir_name="./"
    nums=np.arange(0,2501,1)
    bwv=Asc.Bscan()
    bwv.load(dir_name,nums)
    x1=15.0; dx=-0.5
    tlim=[10,90]
    bwv.set_xaxis(x1,dx,npnt=61)

    fig1=plt.figure(figsize=(6,3.5))
    ax=fig1.add_subplot(111)
    im=bwv.show(ax)
    ax.set_xlim(tlim)
    fsz=12
    ax.tick_params(labelsize=fsz)
    ax.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    ax.set_ylabel("y [mm]",fontsize=fsz)
    plt.colorbar(im)


    #-------------- A-scans -------------
    fig2=plt.figure(figsize=(6,3.5))
    bx=fig2.add_subplot(111)
    bx.grid(True)
    bx.tick_params(labelsize=fsz)
    bwv.normalize()
    awv0=bwv.get_mean()

    num=30
    ycod=bwv.xcod[num]
    awv=bwv.get_ascan(num)

    awv.show_t(bx,name="y="+str(ycod)+"mm")
    awv0.show_t(bx,name="mean",lwd=2)
    bx.set_xlim(tlim)
    bx.legend()
    bx.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    bx.set_ylabel("amplitude",fontsize=fsz)

    #-------------- Group Delay  -------------
    fig3=plt.figure(figsize=(6,4))
    cx=fig3.add_subplot(111)
    cx.grid(True)
    cx.set_xlim([0,2])
    cx.set_ylim([0,30])
    #cx.set_ylim([0,30])
    eps=0.02
    f1=0.2
    f2=1.5
    #awv.get_tdly(ax=cx,eps=eps)
    #awv0.get_tdly(ax=cx,eps=eps)
    awv.show_fft()
    awv0.show_fft()
    awv.show_phi(cx,eps=eps)
    awv0.show_phi(cx,eps=eps)
    #tg=awv.show_dphi(ax=cx,f1=f1,f2=f2,eps=eps)
    #tg0=awv0.show_dphi(ax=cx,f1=f1,f2=f2,eps=eps)
    #print("tg=",tg)
    #print("tg0=",tg0)
    #cx.set_ylim([-4,4])

    plt.show()

    fig1.savefig("bscan.png",bbox_inches="tight")
    fig2.savefig("ascan.png",bbox_inches="tight")
