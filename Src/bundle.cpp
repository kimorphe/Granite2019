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

	char dir_name[128];
	char fname[128],fngrd[128];
	char fnout1[128],fnout2[128],fntmp[128];
	char cbff[128];
	FILE *fp=fopen("bundle.inp","r");

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",dir_name);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fntmp);
	sprintf(fngrd,"%s/%s",dir_name,fntmp);	// Grid file name
	puts(fngrd);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fntmp);
	sprintf(fnout1,"%s/%s",dir_name,fntmp);	// File name of bundled A-scans
	puts(fnout1);

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fntmp);
	sprintf(fnout2,"%s/%s",dir_name,fntmp);	// Output file name (k-vector)
	puts(fnout2);

	Grid Gd;	// Measurement Grid 
	Gd.load(fngrd); // import grid info.
	//---------------------------------------
	sprintf(fname,"%s/scope_%d.csv",dir_name,0);
	Wv1D awv1; // import wave data 1
	awv1.load2(fname,true);	// necessary to get Nt
	int Nf=awv1.FFT(1);

	Array3D WV(Gd.Nx,Gd.Ny,awv1.Nt);	// Measured A-scans
	WV.set_Xa(Gd.Xcod[0],Gd.Ycod[0],awv1.t1);
	WV.set_dx(Gd.dx,Gd.dy,awv1.dt);
	WV.set_Wd();
	WV.print_dim();
	printf(" Loading data from %s\n",dir_name);
	WV.load(dir_name);
	printf(" Writting  data to %s\n",fnout1);
	WV.out(fnout1);


	Array3Dcmplx WVf(Gd.Nx,Gd.Ny,Nf);	// Fourier transform 
	WVf.print_dim();
	WVf.set_Xa(Gd.Xcod[0],Gd.Ycod[0],0.0);
	WVf.set_dx(Gd.dx,Gd.dy,awv1.df);
	WVf.set_Wd();

	printf(" Performing Fourier Transform on A-scans...\n");
	int i,j,k;
	for(i=0;i<Gd.Nx;i++){
		printf("i=%d/%d\n",i,Gd.Nx);
	for(j=0;j<Gd.Ny;j++){
		awv1.amp=WV.A[i][j];
		awv1.mllc=true;
		awv1.Amp=WVf.Z[i][j];
		awv1.FFT(1);
	}
	}
	printf(" Writting Fourier transform to %s\n",fnout2);
	WVf.out(fnout2);
	exit(-1);


	Array3D Phi(Gd.Nx,Gd.Ny,Nf/2);		// Phase(x,y,w)
	char fout[128]="wvf.out";
	WVf.write_zslice(fout,int(1.0/awv1.df));
	for(i=0;i<Gd.Nx;i++){
	for(j=0;j<Gd.Ny;j++){
	for(k=0;k<Nf/2;k++){
		Phi.A[i][j][k]=arg(WVf.Z[i][j][k]);
	}
	}
	}

	Array2D Kx(Gd.Nx,Gd.Ny);		// Kx(x,y,w): x-wave number
	Array2D Ky(Gd.Nx,Gd.Ny);		// Ky(x,y,w): y-wave number
	k=int(1.2/awv1.df);
	double dpx,dpy,wgt;
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);
	wgt=2.0*Gd.dx;
	for(i=0;i<Gd.Nx-1;i++){
	for(j=0;j<Gd.Ny;j++){
		dpx=Phi.A[i+1][j][k]-Phi.A[i][j][k];
		if(dpx> PI) dpx=-(PI2-dpx);
		if(dpx<-PI) dpx=PI2+dpx;
		dpx/=wgt;
		Kx.A[i][j]+=dpx;
		Kx.A[i+1][j]+=dpx;
	}
	}
	for(j=0;j<Gd.Ny;j++){
		Kx.A[0][j]*=2.0;
		Kx.A[Gd.Nx-1][j]*=2.0;
	}

	wgt=2.0*Gd.dy;
	for(i=0;i<Gd.Nx;i++){
	for(j=0;j<Gd.Ny-1;j++){
		dpy=Phi.A[i][j+1][k]-Phi.A[i][j][k];
		if(dpy> PI) dpy=-(PI2-dpy);
		if(dpy<-PI) dpy=PI2+dpy;
		dpy/=wgt;
		Ky.A[i][j]+=dpy;
		Ky.A[i][j+1]+=dpy;
	}
	}
	for(i=0;i<Gd.Nx;i++){
		Ky.A[i][0]*=2.0;
		Ky.A[i][Gd.Ny-1]*=2.0;
	}
	char fnkx[128]="kx.out";
	char fnky[128]="ky.out";
	Kx.out(fnkx);
	Ky.out(fnky);


	exit(-1);
	//WV.Butterworth(0.0,3.0);

	//---------------------------------------
	Array2D Bwv,Tg;
	Bwv=WV.proj();
	puts("done");
	char fn[128]="bwv.out";
	Bwv.out(fn);

	Tg=WV.gdelay(2.9);
	sprintf(fn,"tg.out");
	Tg.out(fn);
	
	return(0);
};
