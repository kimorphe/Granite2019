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

class PHASE{
	public:
		double ***P;	// phase data \phi (x,y,w)
		int ***H;	// histogram h( y,w; dy, dw)
		int Nd[3];
		int Nx,Ny,Nz;	// data
		//	Histogram
		int nf,nt;	// histogram
		double f1,f2;
		double t1,t2;
		int ndat;
		double dx[3];
		double Xa[3],Xb[3],Wd[3];
		void load(char *fn);
		double ***mem_alloc(int nx, int ny, int nz);
		int ***mem_ialloc(int n_y, int n_f, int n_t);
		void dump(char *fn);
		void histogram(double f1, double f2, double t1, double t2, int n_f, int n_t);
		void write_hist0(char *fn);
		void write_hist(char *fn);
	private:
};
void PHASE::write_hist(char *fn){
	FILE *fp=fopen(fn,"w");

	fprintf(fp,"# y1, y2 [mm]\n");
	fprintf(fp,"%lf,%lf\n",Xa[1],Xb[1]);
	fprintf(fp,"# f1, f2 [MHz]\n");
	fprintf(fp,"%lf,%lf\n",f1,f2);
	fprintf(fp,"# t1, t2 [micro sec]\n");
	fprintf(fp,"%lf,%lf\n",t1,t2);
	fprintf(fp,"# Ny, Nf, Nt\n");
	fprintf(fp,"%d,%d,%d\n",Ny,nf,nt);
	fprintf(fp,"# count(y,w,t)\n");
	int i,j,k;
	for(i=0; i<Ny; i++){
	for(j=0; j<nf; j++){
	for(k=0; k<nt; k++){
		fprintf(fp,"%d\n",H[i][j][k]);
	}
	}
	}
	fclose(fp);
};
void PHASE::write_hist0(char *fn){
	FILE *fp=fopen(fn,"w");
	int i,j,k;
	int count;
	double tt;
	double Dt=(t2-t1)/nt;
	double Df=(f2-f1)/nf;
	double Dy=dx[1];
	double ycod,freq,time;
	for(i=0; i<Ny; i++){
		ycod=Xa[1]+Dy*i;
	for(k=0; k<nt; k++){
		time=t1+Dt*k;
		count=0;
	//for(j=0; j<nf; j++){
	for(j=60; j<nf; j++){
		freq=f1+Df*j;
		count+=H[i][j][k];
	}
		fprintf(fp,"%lf %d\n",time,count);
	}
		fprintf(fp,"\n");
	}
	fclose(fp);
};
void PHASE::histogram(double f_1, double f_2, double t_1, double t_2, int n_f, int n_t){

	f1=f_1; f2=f_2; nf=n_f;
	t1=t_1; t2=t_2; nt=n_t;
	H=mem_ialloc(Ny,nf,nt);

	double Dt,Df,Dy;
	Dt=(t2-t1)/(nt-1);
	Df=(f2-f1)/(nf-1);
	Dy=fabs(dx[2]);

	int i,j,k;
	int iw,itof;
	double tof,xcod,ycod,omg,freq;
	double PI2=8.0*atan(1.0);
	for(i=0; i<Nx; i++){
		xcod=Xa[0]+dx[0]*i;
	for(j=0; j<Ny; j++){
		//ycod=Xa[1]+dx[1]*j;
	for(k=0; k<Nz; k++){
		if( P[i][j][k]<0.0) continue;
		freq=Xa[2]+dx[2]*k;
		iw=int((freq-f1)/Df+0.5);
		omg=PI2*freq;
		tof=P[i][j][k]/omg;
		itof=int((tof-t1)/Dt+0.5);
		if(iw <0) continue;
		if(iw >=nf) continue;
		if(itof <0) continue;
		if(itof >=nt) continue;
		H[j][iw][itof]++;
	}
	}
	}
};
void PHASE::dump(char *fn){
	FILE *fp=fopen(fn,"w");
	int i,j,k;
	double freq,omg,ycod;
	double PI2=8.0*atan(1.0);
	for(i=0;i<Nd[0];i++){
	for(j=0;j<Nd[1];j++){
		ycod=Xa[1]+dx[1]*j;
	//for(k=0;k<Nd[2];k++){
	for(k=210;k<220;k++){
		freq=Xa[2]+dx[2]*k;
		omg=PI2*freq;
		if( P[i][j][k]<0.0) continue;
		fprintf(fp,"%lf %lf %lf\n",freq,ycod,P[i][j][k]/omg);
	}
	}
	}
	fclose(fp);
};
void PHASE::load(char *fn){
	char cbff[128];
	FILE *fp=fopen(fn,"r");
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",dx,dx+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Nd,Nd+1);

	for(int i=0;i<3;i++) Xb[i]=Xa[i]+dx[i]*(Nd[i]-1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf, %d\n",Xa+2,dx+2,Nd+2);
	fgets(cbff,128,fp);

	P=PHASE::mem_alloc(Nd[0],Nd[1],Nd[2]);
	int i,j,k;
	double dat[5];
	for(i=0;i<Nd[0];i++){
		printf("i=%d\n",i);
	for(j=0;j<Nd[1];j++){
	for(k=0;k<Nd[2];k++){
		fscanf(fp,"%le, %le, %le, %le, %le \n",dat,dat+1,dat+2,dat+3,dat+4);
		P[i][j][k]=dat[0];
	}
	}
	}
}
double ***PHASE::mem_alloc(int nx, int ny, int nz){
	Nx=nx;
	Ny=ny;
	Nz=nz;

	Nd[0]=Nx;
	Nd[1]=Ny;
	Nd[2]=Nz;
	
	ndat=Nx*Ny*Nz;
	int i,j,k;  

	double *Z3=(double *)malloc(sizeof(double)*ndat);
	for(k=0;k<ndat;k++) Z3[k]=0.0;

	double **Z2=(double **)malloc(sizeof(double *)*nx*ny);	
	for(j=0;j<nx*ny; j++) Z2[j]=Z3+j*nz;

	double ***Z=(double ***)malloc(sizeof(double **)*nx);
	for(i=0;i<nx;i++) Z[i]=Z2+i*ny;

	return(Z);
};
int ***PHASE::mem_ialloc(int nx, int ny, int nz){

	int nsize=nx*ny*nz;
	int i,j,k;  

	int *Z3=(int *)malloc(sizeof(int)*nsize);
	for(k=0;k<nsize;k++) Z3[k]=0;

	int **Z2=(int **)malloc(sizeof(int *)*nx*ny);	
	for(j=0;j<nx*ny; j++) Z2[j]=Z3+j*nz;

	int ***Z=(int ***)malloc(sizeof(int **)*nx);
	for(i=0;i<nx;i++) Z[i]=Z2+i*ny;
	return(Z);

}
//------------
void LinSolve(double **A, double *b, int n){ // Gaussian elimination

        double *y=(double *)malloc(sizeof(double)*n);

        int l,i,j;
        double p,q;

        // Forward Elimination
        for(l=0;l<n; l++){      // row 
                p=A[l][l]; //pivot
                b[l]/=p;
                for(j=l; j<n; j++) A[l][j]/=p;
                for(i=l+1; i<n; i++){
                        q=A[i][l];
                        for(j=l; j<n; ++j) A[i][j]-=A[l][j]*q;
                        b[i]-=b[l]*q;
                }
        } //row

        for(l=n-1;l>=0;l--){
                y[l]=b[l];
                for(j=l+1;j<n;j++) y[l]-=A[l][j]*y[j];
        }

	for(i=0; i<n; i++) b[i]=y[i];
};
void linfit(double x1, double dx, double *y, int n, double *a, double *b){

	int i,j;
	double A[2][2];

	double xb=0.0;
	double yb=0.0;
	double x2b=0.0;
	double xyb=0.0;
	int isum=0;
	double x;
	for(j=0; j<n; j++){
		if(y[j]<0.0) continue;
		isum++;
		x=x1+dx*j;
		xb+=x;
		yb+=y[j];
		x2b+=(x*x);
		xyb+=(x*y[j]);
	}
	xb/=isum;
	yb/=isum;
	xyb/=isum;
	x2b/=isum;

	double Det=x2b-xb*xb;
	*a=(xyb-xb*yb)/Det;
	*b=(-xyb*xb+x2b*yb)/Det;

};

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
void polyfit(double *x, double *y, int n, double *ak, int deg){

	int i,p,q;
	int np=deg+1;
	double **A=mem_alloc_double2d(np,np);
	//double *b=(double *)malloc(sizeof(double)*np);
	for(i=0;i<np;i++) ak[i]=0.0;

	double xp,xq,xpy,xi,yi;

	int count=0;
	for(i=0;i<n;i++){
		xi=x[i];
		yi=y[i];
		if(yi <0.0) continue;
		xp=1.0;
		for(p=0; p<np; p++){
			xq=1.0;
			for(q=0; q<np; q++){
				A[p][q]+=(xp*xq);
				xq*=xi;
			}
			//b[p]+=(xp*yi);
			ak[p]+=(xp*yi);
			xp*=xi;
		}
		count++;
	}


	for(p=0; p<np; p++){
		ak[p]/=count;
		for(q=0; q<np; q++){
			A[p][q]/=count;
		}
	}

	LinSolve(A,ak,np);
	free(A);

};
int main(){
	PHASE PH;
	char fname[128]="phix.out";
	PH.load(fname);
	sprintf(fname,"tofs.dat");
	//PH.dump(fname);
	double f1,f2;
	double t1,t2;
	int nf,nt;
	f1=0.6; f2=2.0; nf=100;
	f1=1.0; f2=1.5; nf=100;
	f1=0.6; f2=1.0; nf=100;
	f1=1.0; f2=1.3; nf=100;
	t1=0.0; t2=10.0; nt=100;
	PH.histogram(f1,f2,t1,t2,nf,nt);

	sprintf(fname,"hist_ywt.dat");
	PH.write_hist(fname);

	return(0);
};
#if DB==0
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
	double omg;
	double PI2=8.0*atan(1.0);
	double C0,C1;
	double *ycod=(double *)malloc(sizeof(double)*Fw.Ny);
	double *Ycod=(double *)malloc(sizeof(double)*Fw.Ny*Fw.Nx);
	double *Zcod=(double *)malloc(sizeof(double)*Fw.Ny*Fw.Nx);
	int deg=2;
	double tofb,tof0,tof1;
	double *ak=(double *)malloc(sizeof(double)*(deg+1));
	for(k=0;k<Fw.Ny;k++) ycod[k]=0.0+Fw.dx[1]*k;
	j=0;
	for(i=0;i<Fw.Nx;i++){
		for(k=0;k<Fw.Ny;k++){
			Ycod[j]=ycod[k];
			j++;
		}
	}
	int l,NN=Fw.Ny*Fw.Nx;
	for(k=nf1;k<=nf2;k++){
		xproj_min(Fws[ksum].Pmin,py[ksum],npy[ksum],0.0,Fw.Nx,Fw.Ny);
		xproj_mean(Fws[ksum].Pmin,pyb[ksum],npyb[ksum],0.0,Fw.Nx,Fw.Ny);
		omg=Fws[ksum].freq*PI2;
		for(j=0;j<Fw.Ny;j++){
			py[ksum][j]/=omg;
			pyb[ksum][j]/=omg;
			fprintf(fp,"%lf %d %lf %d\n",py[ksum][j],npy[ksum][j],pyb[ksum][j],npyb[ksum][j]);
		};

		j=0;
		for(i=0;i<Fw.Nx;i++){
		for(l=0;l<Fw.Ny;l++){
				Zcod[j]=Fws[ksum].Pmin[i][l]/omg;
				j++;
		}
		}
		//fprintf(fp,"\n");
		linfit(0.0,Fw.dx[1],py[ksum],Fw.Ny,&C0,&C1);
		//polyfit(ycod, py[ksum], Fw.Ny, ak, deg);
		polyfit(Ycod, Zcod,NN, ak, deg);
		if(deg==2){
		       	tofb=ak[1]+ak[2]*(Fw.dx[1]*(Fw.Ny-1));
		       	tof0=ak[1];
		       	tof1=ak[1]+2.*ak[2]*(Fw.dx[1]*(Fw.Ny-1));
		}
		if(deg==1) tofb=ak[1];

		printf("%lf %lf %lf %lf %lf\n",omg/PI2,1./C0,1./tofb,1./tof0,1./tof1);
			
		ksum++;
	}
	fclose(fp);
	
	return(0);
}
#endif
