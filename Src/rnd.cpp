#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <complex>
#include <math.h>
#include "kvecs.h"

using namespace std;

int main(){

	std::random_device rd;
	//std::mt19937 mt(rd());
	std::mt19937 mt(11);

	int i,j;

	FSLICE fsl;

	char fname[128];
	int num=92;
	sprintf(fname,"k%d.out",num);
	fsl.load_Grad(fname);
	int Nx=fsl.Nx;
	int Ny=fsl.Ny;
	int ndat=Nx*Ny;
	std::uniform_int_distribution<int> RndI(0,ndat-1);

	int id,ix,iy;
	double kx,ky,kk;
	double xx=0,yy=0;
	double ds=0.5;
	for(j=0;j<100;j++){
			xx=0.0; yy=0.0;
	for(i=0;i<100;i++){
		id=RndI(mt);
		ix=id/Ny;
		iy=id%Ny;
		kx=fsl.Kx[ix][iy];
		ky=fsl.Ky[ix][iy];
		kk=sqrt(kx*kx+ky*ky);
		xx+=(kx/kk*ds);
		yy+=(ky/kk*ds);
		printf("%lf %lf\n",xx,yy);
	};
	printf("\n");
	}


	return(0);
};
