import numpy as np
import matplotlib.pyplot as plt


fname="wvf.out"
fp=open(fname,"r")

dat=fp.readline()
dat=dat.strip().split(",")
Nx=int(dat[0])
Ny=int(dat[1])

print("Nx=",Nx)
print("Ny=",Ny)

Z=[]
for row in fp:
    dat=row.strip().split(",")
    R=float(dat[0])
    I=float(dat[1])
    X=R+1j*I
    Z.append(X)


Z=np.array(Z)
Z=np.reshape(Z,[Nx,Ny])

fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(np.angle(Z),cmap="jet",aspect="equal",interpolation="none",origin="lower")

plt.show()




