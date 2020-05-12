#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "waves.h"

int main(){
	Array2D Kx,Ky;
	char fnamex[128]="kx1.out";
	char fnamey[128]="ky1.out";

	Kx.load(fnamex);
	Ky.load(fnamey);

	//printf("Nx,Ny=%d %d\n",Kx.Nx,Kx.Ny);
	//printf("Nx,Ny=%d %d\n",Ky.Nx,Ky.Ny);

	double xp[2],xi[2];
	int i0,j0,i,j;
	i0=0;
	j0=0;
	int k=0,kmax=400;
	double kx,ky,kk;
	double alph=0.25;
	double ds=Kx.dx[0]*alph;
	int ip,Np=Kx.Nx*2;
	double xp0,yp0;
	xp0=Kx.dx[0]*(i0+0.5);
	yp0=Kx.dx[1]*(j0+0.5);
	for(ip=0;ip<Np;ip++){
		xp[0]=xp0;
		xp[1]=yp0;
		printf("%lf %lf\n",Kx.Xa[0]+xp[0],Kx.Xa[1]+xp[1]);
	for(k=0;k<kmax;k++){
		i=floor(xp[0]/Kx.dx[0]);
		j=floor(xp[1]/Kx.dx[1]);

		if(i<0) break;
		if(j<0) break;
		if(i>=Kx.Nx) break;
		if(j>=Ky.Ny) break;
		kx=-Kx.A[i][j];
		ky=-Ky.A[i][j];
		kk=sqrt(kx*kx+ky*ky);
		xi[0]=kx/kk;
		xi[1]=ky/kk;
		xp[0]+=(xi[0]*ds);
		xp[1]+=(xi[1]*ds);
		printf("%lf %lf\n",Kx.Xa[0]+xp[0],Kx.Xa[1]+xp[1]);
	}
		xp0+=ds;
		printf("\n");
	}
	return(0);
};
