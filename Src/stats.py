import numpy as np
import matplotlib.pyplot as plt

class Stats:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        freq=[]
        kb=[]
        sk=[];
        thb=[]
        sth=[];
        for row in fp:
            dat=row.strip().split(" ")
            freq.append(float(dat[0]))
            kb.append(float(dat[1]))
            sk.append(float(dat[2]))
            thb.append(float(dat[3]))
            sth.append(float(dat[4]))


        self.freq=np.array(freq)
        self.kb=np.array(kb)
        self.thb=np.array(thb)
        self.sk=np.array(sk)
        self.sth=np.array(sth)

        #print(self.thb)

if __name__=="__main__":

    stat=Stats()
    stat.load("mean.out")

    fig1=plt.figure()
    fig2=plt.figure()
    ax=fig1.add_subplot(111)
    bx=fig2.add_subplot(111)
    ax.grid(True)
    bx.grid(True)

    fsz=16
    #ax.set_xlabel("frequency [MHz]",fontsize=fsz)
    #ax.set_ylabel("wave number [mm$^{-1}$]",fontsize=fsz)
    #bx.set_xlabel("frequency [MHz]",fontsize=fsz)
    #bx.set_ylabel("wave direction [deg]",fontsize=fsz)
    ax.tick_params(labelsize=fsz-2)
    bx.tick_params(labelsize=fsz-2)

    lwd=2.0
    stat.thb+=90;

    ax.plot(stat.freq,stat.kb,"b",linewidth=lwd)
    ax.plot(stat.freq,stat.sk,"r",linewidth=lwd)
    bx.plot(stat.freq,stat.thb,"b",linewidth=lwd)
    bx.plot(stat.freq,stat.sth,"r",linewidth=lwd)

    f1=stat.freq[0];
    f2=stat.freq[-1];
    ax.set_xlim([f1,f2])
    bx.set_xlim([f1,f2])

    fig1.savefig("kw_plot.png",bbox_inches="tight")
    fig2.savefig("thw_plot.png",bbox_inches="tight")

    plt.show()
