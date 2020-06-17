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
int main(){
	PHASE PH;
	char cbff[128],fname[128];
	double f1,f2,t1,t2;
	int nf,nt;

	FILE *fp=fopen("tofs.inp","r");

	fgets(cbff,128,fp);
	fscanf(fp,"%s\n",fname);
	printf(" Phase data <-- %s\n",fname);
		PH.load(fname);	// load Fourier phase data
	fgets(cbff,128,fp);
	printf(" Histogram --> %s\n",fname);
	fscanf(fp,"%s\n",fname);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf, %d\n",&f1,&f2,&nf);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf, %d\n",&t1,&t2,&nt);

	fclose(fp);
	//f1=0.6; f2=2.0; nf=100;
	//t1=0.0; t2=10.0; nt=100;

	printf("f1,f2=%lf,%lf [MHz] nf=%d\n",f1,f2,nf);
	printf("t1,t2=%lf,%lf [micro sec] nt=%d\n",t1,t2,nt);
	PH.histogram(f1,f2,t1,t2,nf,nt);

	PH.write_hist(fname);

	return(0);
};
#if DB==0
#endif
