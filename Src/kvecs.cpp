#define DB 2
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "waves.h"

using namespace std;

class FSLICE{
	public:
		double freq;
		double Xa[2],Xb[2],Wd[2],dx[2];
		double **Kx,**Ky,**Phi,**Amp,**E,**F,**G;
		int Nx,Ny,ndat,Nd[2];
		void init(int n, int m);
		void set_Xa(double *xa);
		void set_dx(double *dx);
		void set_Wd();
		void print_domain();
		void Grad();
		void get_slice(complex<double> ***Z,int k);
		void export_Grad(char fn[128]);
		void histogram(double kmin, double kmax, double nbin, double *prob_k, double *prob_a);
	private:
		void malloc_arrays();
		double** mem_alloc(int n, int m);
};
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
	puts("allocating arrays");
	malloc_arrays();
	puts("done !");
};

void FSLICE::malloc_arrays(){
	Kx=mem_alloc(Nx,Ny);
	Ky=mem_alloc(Nx,Ny);
	Phi=mem_alloc(Nx,Ny);
	Amp=mem_alloc(Nx,Ny);
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
	for(i=0;i<Nx-1;i++){ 
	for(j=0;j<Ny;j++){
		dpx=Phi[i+1][j]-Phi[i][j];
		if(dpx> PI) dpx=-(PI2-dpx);
		if(dpx<-PI) dpx=PI2+dpx;
		dpx/=wgt;
		Kx[i][j]+=dpx;
		Kx[i+1][j]+=dpx;
	}
	}

	for(j=0;j<Ny;j++){
		Kx[0][j]*=2.0;
		Kx[Nx-1][j]*=2.0;
	}

	wgt=2.0*dx[1]*PI2;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny-1;j++){
		dpy=Phi[i][j+1]-Phi[i][j];
		if(dpy> PI) dpy=-(PI2-dpy);
		if(dpy<-PI) dpy=PI2+dpy;
		dpy/=wgt;
		Ky[i][j]+=dpy;
		Ky[i][j+1]+=dpy;
	}
	}
	for(i=0;i<Nx;i++){
		Ky[i][0]*=2.0;
		Ky[i][Ny-1]*=2.0;
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
void FSLICE::histogram(double kmin, double kmax, double nbin, double *prob_k, double *prob_a){
	double dk=(kmax-kmin)/nbin;
	double da=360./nbin;
	double PI=4.0*atan(1.0);
	double xi,xi0,xi1,alph;
	double count=1.0/ndat;
	double asum,amp;
	int i,j,ibin;
	double k_mean,k_sig,a_mean,a_sig;
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
	//fprintf(fout,"%lf %lf %lf %lf %lf\n",freq,k_mean,k_sig,a_mean,a_sig);
};

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
	double tmp=1.5;
	printf("freq=%lf, indx=%d\n",tmp,WVf.get_index(tmp,2));
	// -------------------------------
	
	FSLICE Fw;

	Fw.init(WVf.Nx,WVf.Ny);
	Fw.set_Xa(WVf.Xa);
	Fw.set_dx(WVf.dx);
	Fw.set_Wd();
	Fw.print_domain();

	char fnkx[128],fnky[128];
	int i,j,k;
	int k1,k2,inc,ksum;
	inc=int(df/WVf.dx[2]);
	if(inc<1) inc=1;
	k1=WVf.get_index(f1,2);
	k2=WVf.get_index(f2,2);

	ksum=(k2-k1)/inc+1;
	printf("ksum=%d\n",ksum);
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

	Array2D prob_k(ksum,nbin);
	Array2D prob_a(ksum,nbin);

	Prob_k.set_dx(dk,df);
	Prob_a.set_dx(da,df);
	Prob_k.set_Xa(kmin,f1);
	Prob_a.set_Xa(-270.0,f1);

	ksum=0;
	for(k=k1;k<=k2;k+=inc){
		Fw.freq=WVf.dx[2]*k+WVf.Xa[2];
		printf("%d f=%lf, %lf[MHz]\n",ksum,Fw.freq,k*WVf.dx[2]);
		Fw.get_slice(WVf.Z,k);	
		sprintf(fnkx,"k%d.out",ksum);
		Fw.Grad();
		Fw.export_Grad(fnkx);
		Fw.histogram(kmin,kmax,nbin,prob_k.A[ksum],prob_a.A[ksum]);
		ksum++;
	};

	Array3D Phi(WVf.Nx,WVf.Ny,WVf.Nz); // Phase(x,y,w)
	Array2D Kx(WVf.Nx, WVf.Ny);	// Kx(x,y,w): x-wave number
	Array2D Ky(WVf.Nx, WVf.Ny);	// Ky(x,y,w): y-wave number

	Kx.set_Xa(WVf.Xa[0],WVf.Xa[1]);
	Kx.set_dx(WVf.dx[0],WVf.dx[1]);
	Kx.set_Wd();
	Ky.set_Xa(WVf.Xa[0],WVf.Xa[1]);
	Ky.set_dx(WVf.dx[0],WVf.dx[1]);
	Ky.set_Wd();

	int ndat=WVf.Nx*WVf.Ny;

	for(i=0;i<WVf.Nx;i++){
	for(j=0;j<WVf.Ny;j++){
	for(k=0;k<WVf.Nz;k++){
		Phi.A[i][j][k]=arg(WVf.Z[i][j][k]);
	}
	}
	}

	if(df < WVf.dx[2]) df=WVf.dx[2];
	double dpx,dpy,wgt;
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);

	//inc=int(df/WVf.dx[2]);
	//k1=int(f1/WVf.dx[2]);
	//k2=int(f2/WVf.dx[2]);
	double freq=f1;
	//while(freq<=f2){

	f1=k1*WVf.dx[2];
	f2=k2*WVf.dx[2];
	df=inc*WVf.dx[2];



	FILE *fout=fopen("mean.out","w");
	fprintf(fout,"# freq.[MHz], <k>, sig_k, <th>, sig_th\n");
	double k_mean,k_sig,a_mean,a_sig;
	ksum=0;
	double count=1.0/ndat;
	for(k=k1;k<=k2;k+=inc){
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
