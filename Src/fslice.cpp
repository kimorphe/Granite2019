#define DB 2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "waves.h"
#include "kvecs.h"

using namespace std;
double** FSLICE::mem_alloc(int n, int m){
	int i;
	double **A,*A2;
	A2=(double *)malloc(sizeof(double)*n*m);
	A=(double **)malloc(sizeof(double*)*n);
	for(i=0;i<n*m;i++) A2[i]=0.0;
	for(i=0;i<n;i++)  A[i]=A2+m*i;
	return(A);
};
void FSLICE::init(int n, int m){
	Nx=n; Ny=m; ndat=n*m;
	Nd[0]=n; Nd[1]=m;
	malloc_arrays();
};

void FSLICE::malloc_arrays(){
	Kx=mem_alloc(Nx,Ny);
	Ky=mem_alloc(Nx,Ny);
	Kxx=mem_alloc(Nx,Ny);
	Kyy=mem_alloc(Nx,Ny);
	Kxy=mem_alloc(Nx,Ny);
	Phi=mem_alloc(Nx,Ny);
	Amp=mem_alloc(Nx,Ny);
	is_alloc=1;
};
void FSLICE::set_Xa(double *xa){
	Xa[0]=xa[0];
	Xa[1]=xa[1];
};
void FSLICE::set_dx(double *h){
	dx[0]=h[0];
	dx[1]=h[1];
};
void FSLICE::set_Wd(){
	for(int i=0;i<2;i++){
		Wd[i]=dx[i]*(Nd[i]-1);
		Xb[i]=Xa[i]+(Wd[i]-1);
	}
};
void FSLICE::print_domain(){
	puts("------      DOMAIN SIZE ------- ");
	printf("Xa=%lf %lf\n",Xa[0],Xa[1]);
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);
	printf("Nd=%d %d\n",Nd[0],Nd[1]);
	printf("dx=%lf %lf\n",dx[0],dx[1]);
	puts("------------------------------- ");
};
void FSLICE::get_slice(complex<double> ***Z,int k){
	int i,j;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		Phi[i][j]=arg(Z[i][j][k]);
		Amp[i][j]=abs(Z[i][j][k]);
	}
	}
};

void FSLICE::Grad(){
		
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);
	double wgt=2.0*dx[0]*PI2;
	double dpx,dpy;
	int i,j;
	for(i=0;i<Nx;i++){ 
	for(j=0;j<Ny;j++){
		Kx[i][j]=0.0; 
		Ky[i][j]=0.0; 
		Kxx[i][j]=0.0;
		Kyy[i][j]=0.0;
		Kxy[i][j]=0.0;
	}
	}

	double dx2=dx[0]*dx[0];
	double dy2=dx[1]*dx[1];

	for(i=0;i<Nx-1;i++){ 
	for(j=0;j<Ny;j++){
		dpx=Phi[i+1][j]-Phi[i][j];
		if(dpx> PI) dpx=-(PI2-dpx);
		if(dpx<-PI) dpx=PI2+dpx;
		Kxx[i][j]+=(dpx/dx2);
		Kxx[i+1][j]-=(dpx/dx2);

		dpx/=wgt;
		Kx[i][j]+=dpx;
		Kx[i+1][j]+=dpx;

		//Kyx[i][j]+=dpx
	}
	}

	for(j=0;j<Ny;j++){
		Kx[0][j]*=2.0;
		Kx[Nx-1][j]*=2.0;

		Kxx[0][j]=0.0;
		Kxx[Nx-1][j]=0.0;
	}

	wgt=2.0*dx[1]*PI2;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny-1;j++){
		dpy=Phi[i][j+1]-Phi[i][j];
		if(dpy> PI) dpy=-(PI2-dpy);
		if(dpy<-PI) dpy=PI2+dpy;
		Kyy[i][j]+=(dpy/dy2);
		Kyy[i][j+1]-=(dpy/dy2);

		dpy/=wgt;
		Ky[i][j]+=dpy;
		Ky[i][j+1]+=dpy;
	}
	}
	for(i=0;i<Nx;i++){
		Ky[i][0]*=2.0;
		Ky[i][Ny-1]*=2.0;

		Kyy[i][0]=0.0;
		Kyy[i][Ny-1]=0.0;
	}

	for(i=1;i<Nx-1;i++){
	for(j=1;j<Ny-1;j++){
		Kxy[i][j]=(Kx[i][j+1]-Kx[i][j-1])/(2.*dx[1]);
		Kxy[i][j]=Kxy[i][j]+(Ky[i+1][j]-Ky[i-1][j])/(2.*dx[0]);
		Kxy[i][j]*=0.5;
	}
	}


};
void FSLICE::export_Grad(char fn[128]){
	FILE *fp=fopen(fn,"w");

	int i,j;
	fprintf(fp,"# frequency [MHz]\n");
	fprintf(fp,"%lf\n",freq);
	fprintf(fp,"# Nx, Ny\n");
	fprintf(fp,"%d,%d\n",Nx,Ny);
	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf,%lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# dx[0:1]\n");
	fprintf(fp,"%lf,%lf\n",dx[0],dx[1]);
	fprintf(fp,"# kx, ky (for x{ for y})\n");
	for( i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fprintf(fp,"%le,%le\n",Kx[i][j],Ky[i][j]);
	}
	}
	fclose(fp);
};
void FSLICE::load_Grad(char fn[128]){
	FILE *fp=fopen(fn,"r");
	int i,j;
	char cbff[128];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&freq);
	printf("#freq=%lf\n",freq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d,%d\n",&Nx,&Ny);
	printf("#Nx,Ny=%d,%d\n",Nx,Ny);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",dx,dx+1);
	fprintf(fp,"# kx, ky (for x{ for y})\n");
	fgets(cbff,128,fp);
	double kx,ky;
	if(is_alloc==0){
		Kx=mem_alloc(Nx,Ny);
		Ky=mem_alloc(Nx,Ny);
	};
	for( i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fscanf(fp,"%le,%le\n",&kx, &ky);
		Kx[i][j]=kx;
		Ky[i][j]=ky;
	}
	}
	fclose(fp);
};
void FSLICE::export_Hess(char fn[128]){
	FILE *fp=fopen(fn,"w");
	int i,j;
	fprintf(fp,"# frequency [MHz]\n");
	fprintf(fp,"%lf\n",freq);
	fprintf(fp,"# Nx, Ny\n");
	fprintf(fp,"%d,%d\n",Nx,Ny);
	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf,%lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# dx[0:1]\n");
	fprintf(fp,"%lf,%lf\n",dx[0],dx[1]);
	fprintf(fp,"# kx, ky (for x{ for y})\n");
	for( i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fprintf(fp,"%le,%le,%le\n",Kxx[i][j],Kyy[i][j],Kxy[i][j]);
	}
	}
	fclose(fp);
};
void FSLICE::histogram(double kmin, double kmax, double nbin, double *prob_k, double *prob_a){
	double dk=(kmax-kmin)/nbin;
	double da=360./nbin;
	double PI=4.0*atan(1.0);
	double xi,xi0,xi1,alph;
	double count=1.0/ndat;
	double asum,amp;
	int i,j,ibin;
	k_mean=0.0; 
	k_sig=0.0; 
	a_mean=0.0; 
	a_sig=0.0;
	asum=0.0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		xi0=Kx[i][j];
		xi1=Ky[i][j];
		xi=sqrt(xi0*xi0+xi1*xi1);
		alph=asin(xi1/xi);
		if(xi0<0.0) alph=-PI-alph;
		alph=alph/PI*180.0;

		ibin=(xi-kmin)/dk;
		//if(ibin>=0 && ibin <nbin) Prob_k.A[ibin][ksum]+=count; 
		if(ibin>=0 && ibin <nbin) prob_k[ibin]+=count; 

		//ibin=(alph+180.0)/da;
		ibin=(alph+270.0)/da;
		if(ibin>=0 && ibin <nbin) prob_a[ibin]+=count;

		amp=Amp[i][j];
		amp*=amp;
		asum+=amp;

		k_mean+=(xi*amp);
		k_sig+=(xi*xi*amp);
		a_mean+=(alph*amp);
		a_sig+=(alph*alph*amp);
	}
	}

	k_mean/=asum;
	k_sig/=asum;
	a_mean/=asum;
	a_sig/=asum;
	k_sig=sqrt(k_sig-k_mean*k_mean);
	a_sig=sqrt(a_sig-a_mean*a_mean);
};

