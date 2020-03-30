#define DB 2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
//#include "fft.h"
//
#include "waves.h"

using namespace std;

int main(){

	Array3Dcmplx WVf;

	char cbff[128];
	char dir_name[128],fname[128],fntmp[128];
	FILE *fp=fopen("kvecs.inp","r");

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",dir_name);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fntmp);
	fgets(cbff,128,fp);
	double f1,f2,df;
	fscanf(fp,"%lf, %lf, %lf\n",&f1,&f2,&df);
	printf("(f1,f2,df=%lf, %lf, %lf\n",f1,f2,df);
	fclose(fp);

	sprintf(fname,"%s/%s",dir_name,fntmp);
	printf("Loading data from %s\n",fname);
	WVf.load(fname);



	//char fnout1[128],fnout2[128];

	Array3D Phi(WVf.Nx,WVf.Ny,WVf.Nz); // Phase(x,y,w)
	Array2D Kx(WVf.Nx, WVf.Ny);	// Kx(x,y,w): x-wave number
	Array2D Ky(WVf.Nx, WVf.Ny);	// Ky(x,y,w): y-wave number

	int i,j,k;
	for(i=0;i<WVf.Nx;i++){
	for(j=0;j<WVf.Ny;j++){
	for(k=0;k<WVf.Nz;k++){
		Phi.A[i][j][k]=arg(WVf.Z[i][j][k]);
	}
	}
	}

	if(df < WVf.dx[2]) df=WVf.dx[2];
	int k1,k2,inc,ksum;
	double dpx,dpy,wgt;
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);
	char fnkx[128],fnky[128];

	inc=int(df/WVf.dx[2]);
	k1=int(f1/WVf.dx[2]);
	k2=int(f2/WVf.dx[2]);
	ksum=0;
	double freq=f1;
	while(freq<=f2){
		k=int(freq/WVf.dx[2]);
		printf("f=%lf, %lf[MHz]\n",freq,k*WVf.dx[2]);
		Kx.freq=k*WVf.dx[2];
		Ky.freq=k*WVf.dx[2];
		wgt=2.0*WVf.dx[0];
		for(i=0;i<WVf.Nx-1;i++){ 
		for(j=0;j<WVf.Ny;j++){
			dpx=Phi.A[i+1][j][k]-Phi.A[i][j][k];
			if(dpx> PI) dpx=-(PI2-dpx);
			if(dpx<-PI) dpx=PI2+dpx;
			dpx/=wgt;
			Kx.A[i][j]+=dpx;
			Kx.A[i+1][j]+=dpx;
		}
		}
		for(j=0;j<WVf.Ny;j++){
			Kx.A[0][j]*=2.0;
			Kx.A[WVf.Nx-1][j]*=2.0;
		}

		wgt=2.0*WVf.dx[1];
		for(i=0;i<WVf.Nx;i++){
		for(j=0;j<WVf.Ny-1;j++){
			dpy=Phi.A[i][j+1][k]-Phi.A[i][j][k];
			if(dpy> PI) dpy=-(PI2-dpy);
			if(dpy<-PI) dpy=PI2+dpy;
			dpy/=wgt;
			Ky.A[i][j]+=dpy;
			Ky.A[i][j+1]+=dpy;
		}
		}
		for(i=0;i<WVf.Nx;i++){
			Ky.A[i][0]*=2.0;
			Ky.A[i][WVf.Ny-1]*=2.0;
		}
		sprintf(fnkx,"kx%d.out",ksum);
		sprintf(fnky,"ky%d.out",ksum);
		Kx.out(fnkx);
		Ky.out(fnky);
		Kx.clear();
		Ky.clear();
		ksum++;
		freq+=df;
	}
	//WV.Butterworth(0.0,3.0);
	
	return(0);
}
