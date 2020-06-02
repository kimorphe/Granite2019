#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <complex>
#include "kvecs.h"

using namespace std;

int main(){
	int ndat=100;

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> RndI(0,ndat-1);

	int i;
	//for(i=0;i<ndat;i++){
//		printf("%d %d\n",i,RndI(mt));
	//}

	FSLICE fsl;

	char fname[128];
	int num=80;
	sprintf(fname,"k%d.out",num);
	fsl.load_Grad(fname);
	int Nx=fsl.Nx;
	int Ny=fsl.Ny;
	printf("%d %d, %lf %lf\n",Nx,Ny,fsl.Kx[Nx-1][Ny-1],fsl.Ky[Nx-1][Ny-1]);

	return(0);
};
