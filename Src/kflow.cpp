#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "waves.h"

int main(){
	Array2D Kx,Ky;
	char fnamex[128]="kx100.out";
	char fnamey[128]="ky100.out";

	Kx.load(fnamex);
	Ky.load(fnamey);

	//printf("Nx,Ny=%d %d\n",Kx.Nx,Kx.Ny);
	//printf("Nx,Ny=%d %d\n",Ky.Nx,Ky.Ny);

	double xp[2],xi[2];
	int i0,j0,i,j;
	i0=0;
	j0=0;
	int k=0,kmax=1000;
	double kx[4],ky[4],phi[4],kk;
	double et[2];
	double alph=0.20;
	double ds=Kx.dx[0]*alph;
	int ip,Np=Kx.Nx*5;
	double xp0,yp0;
	double dx[2],Xa[2];

	Xa[0]=Kx.Xa[0]; Xa[1]=Kx.Xa[1];
	dx[0]=Kx.dx[0]; dx[1]=Kx.dx[1];
	xp0=dx[0]*(i0+0.5)+Xa[0];
	yp0=dx[1]*(j0+0.5)+Xa[1];
	int l;
	for(ip=0;ip<Np;ip++){
		xp[0]=xp0;
		xp[1]=yp0;
		printf("%lf %lf\n",xp[0],xp[1]);
	for(k=0;k<kmax;k++){
		i=floor((xp[0]-Xa[0])/dx[0]);
		j=floor((xp[1]-Xa[1])/dx[1]);

		if(i<0) break;
		if(j<0) break;
		if(i>=Kx.Nx-1) break;
		if(j>=Ky.Ny-1) break;
		kx[0]=Kx.A[i][j];
		kx[1]=Kx.A[i+1][j];
		kx[2]=Kx.A[i+1][j+1];
		kx[3]=Kx.A[i][j+1];

		ky[0]=Ky.A[i][j];
		ky[1]=Ky.A[i+1][j];
		ky[2]=Ky.A[i+1][j+1];
		ky[3]=Ky.A[i][j+1];

		et[0]=2.*((xp[0]-Xa[0])/dx[0]-i)-1.0;
		et[1]=2.*((xp[1]-Xa[1])/dx[1]-j)-1.0;
		phi[0]=0.25*(1.0-et[0])*(1.0-et[1]);
		phi[1]=0.25*(1.0+et[0])*(1.0-et[1]);
		phi[2]=0.25*(1.0+et[0])*(1.0+et[1]);
		phi[3]=0.25*(1.0-et[0])*(1.0+et[1]);
		xi[0]=0.0;
		xi[1]=0.0;
		//phi[0]=1.0;phi[1]=0.0;phi[2]=0.0;phi[3]=0.0;
		for(l=0;l<4;l++){
			xi[0]+=(phi[l]*kx[l]);
			xi[1]+=(phi[l]*ky[l]);
		}
		kk=sqrt(xi[0]*xi[0]+xi[1]*xi[1]);
		xi[0]/=kk;
		xi[1]/=kk;
		xp[0]+=(xi[0]*abs(ds));
		xp[1]+=(xi[1]*abs(ds));
		printf("%lf %lf\n",xp[0],xp[1]);
	}
		xp0+=ds;
		printf("#\n");
	}
	return(0);
};
