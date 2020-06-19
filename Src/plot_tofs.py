#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt


class tStats:
    def load(self,fname):
        fp=open(fname,"r")
        cbff=fp.readline()
        print(cbff)
        ycod=[]
        tsig=[]
        tave=[]
        tmax=[]
        for row in fp:
            dat=row.strip().split(",")
            ycod.append(float(dat[0]))
            tsig.append(float(dat[1]))
            tave.append(float(dat[2]))
            tmax.append(float(dat[3]))

        self.ycod=np.array(ycod)
        self.tsig=np.array(tsig)
        self.tave=np.array(tave)
        self.tmax=np.array(tmax)

        fp.close()

if __name__=="__main__":
    fname="tstat.dat"

    Tf=tStats()
    Tf.load(fname)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.grid(True)

    ax.plot(Tf.tmax,-Tf.ycod)
    ax.plot(Tf.tave,-Tf.ycod)
    ax.plot(Tf.tave+Tf.tsig,-Tf.ycod,"--")
    ax.plot(Tf.tave-Tf.tsig,-Tf.ycod,"--")

    fsz=14
    ax.tick_params(labelsize=fsz)

    plt.show()


