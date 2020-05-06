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

	Array3D Phi(WVf.Nx,WVf.Ny,WVf.Nz); // Phase(x,y,w)
	Array2D Kx(WVf.Nx, WVf.Ny);	// Kx(x,y,w): x-wave number
	Array2D Ky(WVf.Nx, WVf.Ny);	// Ky(x,y,w): y-wave number

	int ndat=WVf.Nx*WVf.Ny;

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
	if(inc<1) inc=1;
	double freq=f1;
	//while(freq<=f2){
	ksum=(k2-k1)/inc+1;
	printf("ksum=%d\n",ksum);

	f1=k1*WVf.dx[2];
	f2=k2*WVf.dx[2];
	df=inc*WVf.dx[2];


	double kmin,kmax,dk;
	int nbin,ibin;
	double xi,xi0,xi1,alph,da;

	kmin=0;
	kmax=1.2;
	nbin=50;
	dk=(kmax-kmin)/nbin;
	da=360./nbin;
	Array2D Prob_k(nbin,ksum);	// |K|(x,y,w)
	Array2D Prob_a(nbin,ksum);	// |K|(x,y,w)

	Prob_k.set_dx(dk,df);
	Prob_a.set_dx(da,df);
	Prob_k.set_Xa(kmin,f1);
	//Prob_a.set_Xa(-180.0,f1);
	Prob_a.set_Xa(-270.0,f1);

	FILE *fout=fopen("mean.out","w");
	fprintf(fout,"# freq.[MHz], <k>, sig_k, <th>, sig_th\n");
	double k_mean,k_sig,a_mean,a_sig;
	ksum=0;
	double count=1.0/ndat;
	for(k=k1;k<=k2;k+=inc){
		//k=int(freq/WVf.dx[2]);
		printf("%d f=%lf, %lf[MHz]\n",ksum,freq,k*WVf.dx[2]);
		Kx.freq=k*WVf.dx[2];
		Ky.freq=k*WVf.dx[2];
		wgt=2.0*WVf.dx[0]*PI2;
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

		wgt=2.0*WVf.dx[1]*PI2;
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
		
		k_mean=0.0; k_sig=0.0;
		a_mean=0.0; a_sig=0.0;
		double asum,Amp;
		asum=0.0;
		for(i=0;i<WVf.Nx;i++){
		for(j=0;j<WVf.Ny;j++){
			xi0=Kx.A[i][j];
			xi1=Ky.A[i][j];
			xi=sqrt(xi0*xi0+xi1*xi1);
			//alph=acos(xi0/xi);
			//if(xi1<0.0) alph=-alph;
			//alph=alph/PI*180.0;
			alph=asin(xi1/xi);
			if(xi0<0.0) alph=-PI-alph;
			alph=alph/PI*180.0;

			ibin=(xi-kmin)/dk;
			if(ibin>=0 && ibin <nbin) Prob_k.A[ibin][ksum]+=count; 

			//ibin=(alph+180.0)/da;
			ibin=(alph+270.0)/da;
			if(ibin>=0 && ibin <nbin) Prob_a.A[ibin][ksum]+=count;

			Amp=abs(WVf.Z[i][j][k]);
			Amp*=Amp;
			asum+=Amp;

			k_mean+=(xi*Amp);
			k_sig+=(xi*xi*Amp);
			a_mean+=(alph*Amp);
			a_sig+=(alph*alph*Amp);
		}
		}

		k_mean/=asum;
		k_sig/=asum;
		a_mean/=asum;
		a_sig/=asum;

		k_sig=sqrt(k_sig-k_mean*k_mean);
		a_sig=sqrt(a_sig-a_mean*a_mean);
		fprintf(fout,"%lf %lf %lf %lf %lf\n",freq,k_mean,k_sig,a_mean,a_sig);

		Kx.clear();
		Ky.clear();
		ksum++;
		freq+=df;
	}
	//WV.Butterworth(0.0,3.0);
	sprintf(fnkx,"Klen.out");
	Prob_k.out(fnkx);
	sprintf(fnky,"Kalp.out");
	Prob_a.out(fnky);
	
	return(0);
}
