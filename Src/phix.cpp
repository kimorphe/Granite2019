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
int main(){

	Array3Dcmplx WVf;

	char cbff[128];
	char dir_name[128],fname[128],fntmp[128];


	FILE *fp=fopen("kvecs.inp","r");

//   Read General Input Data
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
	int Nf=WVf.Nz;

// -------------------------------
	
	FSLICE Fw, *Fws;
	Fws=(FSLICE *)malloc(sizeof(FSLICE)*Nf);

	int i,j,k;

	// Histogram Parameters 
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
	//for(k=k1;k<=k2;k+=inc){	
	printf("Nf=%d\n",Nf);
	int ksum=0;
	//for(k=0; k<Nf; k++){	
	//for(k=280; k<281; k++){	
	for(k=180; k<181; k++){	
		Fw.init(WVf.Nx,WVf.Ny);
		Fw.set_Xa(WVf.Xa);
		Fw.set_dx(WVf.dx);
		Fw.set_Wd();

		Fw.freq=WVf.dx[2]*k+WVf.Xa[2];	// set frequency [MHz]
		printf("%lf\n",Fw.freq);
		Fw.get_slice(WVf.Z,k);	// get Fourier transform
		Fw.Integrate2();
		//Fw.Grad();	// evaluate k-vector field
		//Fw.histogram(kmin,kmax,nbin,prob_k.A[ksum],prob_a.A[ksum]); // get histogram
		//fprintf(fout,"%lf %lf %lf %lf %lf\n",Fw.freq,Fw.k_mean,Fw.k_sig,Fw.a_mean,Fw.a_sig);
		Fws[k]=Fw;
		ksum++;
	};
	exit(-1);
	char fnkx[128],fnky[128];
	sprintf(fnkx,"Klen.out");
	sprintf(fnky,"Kalp.out");
	prob_k.out(fnkx);
	prob_a.out(fnky);
	puts(fnky);

	char fnout[128]="kvec.out";
	fp=fopen(fnout,"w");
	fprintf(fp,"# Xa[0:2]\n");
	fprintf(fp,"%lf, %lf\n",Fw.Xa[0],Fw.Xa[1]);
	fprintf(fp,"# dx[0:2]\n");
	fprintf(fp,"%lf, %lf\n",Fw.dx[0],Fw.dx[1]);
	fprintf(fp,"# Nd[0:2]\n");
	fprintf(fp,"%d, %d\n",Fw.Nx,Fw.Ny);
	fprintf(fp,"# f1, df, Nf\n");
	fprintf(fp,"%lf, %lf, %d\n",WVf.Xa[2],WVf.dx[2],WVf.Nd[2]);
	fprintf(fp,"# wave numb.(kx,ky)\n");
	for(i=0;i<Fw.Nx;i++){
	for(j=0;j<Fw.Ny;j++){
	for(k=0;k<Nf;k++){
		fprintf(fp,"%lf,%f\n",Fws[k].Kx[i][j],Fws[k].Ky[i][j]);
	}
	}
	}
	fclose(fp);
	printf("nline=%d\n",Fw.Nx*Fw.Ny*ksum);
	
	
	return(0);
}
