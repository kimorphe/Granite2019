#define DB 2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "waves.h"
#include "kvecs.h"
#include "heap.h"

using namespace std;

//---------------------------------------------------------------
double **mem_alloc_double2d(int nx, int ny){
	int ndat=nx*ny;
	double *p=(double *)malloc(sizeof(double)*ndat);
	double **A=(double **)malloc(sizeof(double*)*nx);
	int i;
	for(i=0;i<ndat;i++) p[i]=0.0;
	for(i=0;i<nx;i++) A[i]=p+i*ny;
	return(A);
};
int **mem_alloc_int2d(int nx, int ny){
	int ndat=nx*ny;
	int *p=(int *)malloc(sizeof(int)*ndat);
	int **A=(int **)malloc(sizeof(int*)*nx);
	int i;
	for(i=0;i<ndat;i++) p[i]=0;
	for(i=0;i<nx;i++) A[i]=p+i*ny;
	return(A);
}
int main(){

	Array3Dcmplx WVf;

	char cbff[128];
	char dir_name[128],fname[128],fntmp[128];


	FILE *fp=fopen("phix.inp","r");

//   	Read General Input Data
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",dir_name);
	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fntmp);
	fgets(cbff,128,fp);
	double f1,f2,df;
	fscanf(fp,"%lf, %lf\n",&f1,&f2);
	printf("(f1,f2)=%lf,%lf\n",f1,f2);
	fclose(fp);

	sprintf(fname,"%s/%s",dir_name,fntmp);
	printf("Loading data from %s\n",fname);
	WVf.load(fname);
	int Nf=WVf.Nz;

// -------------------------------
	
	FSLICE Fw,*Fws;
	Fws=(FSLICE *)malloc(sizeof(FSLICE)*Nf);

	int nf1=WVf.get_index(f1,2);
	int nf2=WVf.get_index(f2,2);
	f1=WVf.get_cod(nf1,2);
	f2=WVf.get_cod(nf2,2);
	printf("f1,f2=%lf,%lf\n",f1,f2);

	int i,j,k;

	// Histogram Parameters 
/*
	int nbin,ibin;
	double kmin,kmax,dk;
	double xi,xi0,xi1,alph,da;
	kmin=0; 
	kmax=1.2;
	nbin=50;
	dk=(kmax-kmin)/nbin;
	da=360./nbin;
	Array2D prob_k(Nf,nbin);
	Array2D prob_a(Nf,nbin);
	prob_k.set_dx(df,dk);
	prob_a.set_dx(df,da);
	prob_k.set_Xa(f1,kmin);
	prob_a.set_Xa(f1,-270.0);
	FILE *fout=fopen("mean.out","w");
	fprintf(fout,"# freq.[MHz], <k>, sig_k, <th>, sig_th\n");
*/
	int ksum=0;
	//for(k=380; k<381; k++){	
	for(k=nf1;k<=nf2;k++){
		Fw.init(WVf.Nx,WVf.Ny);
		Fw.set_Xa(WVf.Xa);
		Fw.set_dx(WVf.dx);
		Fw.set_Wd();

		//Fw.freq=WVf.dx[2]*k+WVf.Xa[2];	// set frequency [MHz]
		Fw.freq=WVf.get_cod(k,2);
		printf("frq=%lf[MHz]\n",Fw.freq);
		Fw.get_slice(WVf.Z,k);	// get Fourier transform
		Fw.Integrate2();
		//Fw.Grad();	// evaluate k-vector field
		//Fw.histogram(kmin,kmax,nbin,prob_k.A[ksum],prob_a.A[ksum]); // get histogram
		//fprintf(fout,"%lf %lf %lf %lf %lf\n",Fw.freq,Fw.k_mean,Fw.k_sig,Fw.a_mean,Fw.a_sig);
		Fws[ksum]=Fw;
		ksum++;
	};
	/*
	char fnkx[128],fnky[128];
	sprintf(fnkx,"Klen.out");
	sprintf(fnky,"Kalp.out");
	prob_k.out(fnkx);
	prob_a.out(fnky);
	puts(fnky);
	*/

	char fnout[128]="phix.out";
	fp=fopen(fnout,"w");
	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf, %lf\n",Fw.Xa[0],Fw.Xa[1]);
	fprintf(fp,"# dx[0:1]\n");
	fprintf(fp,"%lf, %lf\n",Fw.dx[0],Fw.dx[1]);
	fprintf(fp,"# Nd[0:1]\n");
	fprintf(fp,"%d, %d\n",Fw.Nx,Fw.Ny);
	fprintf(fp,"# f1, df, Nf\n");
	//fprintf(fp,"%lf, %lf, %d\n",WVf.Xa[2],WVf.dx[2],WVf.Nd[2]);
	fprintf(fp,"%lf, %lf, %d\n",f1,WVf.dx[2],ksum);
	fprintf(fp,"# phase (min, max, mean, var) \n");
	double pmin,pave,pmax,pvar,ppnt;
	for(i=0;i<Fw.Nx;i++){
	for(j=0;j<Fw.Ny;j++){
		ksum=0;
	//for(k=0;k<Nf;k++){
	for(k=nf1;k<=nf2;k++){
		pmin=Fws[ksum].Pmin[i][j];
		pmax=Fws[ksum].Pmax[i][j];
		pave=Fws[ksum].Pave[i][j];
		pvar=Fws[ksum].Pvar[i][j];
		ppnt=Fws[ksum].Ppnt[i][j];

		fprintf(fp,"%lf,%f,%lf,%lf,%lf\n",pmin,pmax,pave,pvar,ppnt);
		ksum++;
	}
	}
	}
	fclose(fp);
	//printf("nline=%d\n",Fw.Nx*Fw.Ny*ksum);
	//for(k=0;k<Nf;k++){
	ksum=0;
	int nf12=nf2-nf1+1;
	int **npy=mem_alloc_int2d(nf12,Fw.Ny);
	double **py=mem_alloc_double2d(nf12,Fw.Ny);
	int **npyb=mem_alloc_int2d(nf12,Fw.Ny);
	double **pyb=mem_alloc_double2d(nf12,Fw.Ny);

	fp=fopen("tofs.out","w");
	//fprintf(fp,"# Nf, Ny\n");
	//fprintf(fp,"%d, %d\n",nf12,Fw.Ny);
	fprintf(fp,"# f1, df, Nf\n");
	fprintf(fp,"%lf, %lf, %d\n",f1,WVf.dx[2],nf12);
	fprintf(fp,"# y1, dy, Ny\n");
	fprintf(fp,"%lf, %lf, %d\n",Fw.Xa[1],Fw.dx[1],Fw.Ny);
	fprintf(fp,"# phase(min, count_min, mean, count_mean)\n");
	for(k=nf1;k<=nf2;k++){
		xproj_min(Fws[ksum].Pmin,py[ksum],npy[ksum],0.0,Fw.Nx,Fw.Ny);
		xproj_mean(Fws[ksum].Pmin,pyb[ksum],npyb[ksum],0.0,Fw.Nx,Fw.Ny);

		for(j=0;j<Fw.Ny;j++){
			fprintf(fp,"%lf %d %lf %d\n",py[ksum][j],npy[ksum][j],pyb[ksum][j],npyb[ksum][j]);
		};
		//fprintf(fp,"\n");
		ksum++;

	}
	fclose(fp);
	
	return(0);
}
