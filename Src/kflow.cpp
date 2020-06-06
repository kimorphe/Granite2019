#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <random>
#include "fft.h"

//#include "waves.h"
//#include "kvecs.h"
void FFT2D( complex<double> **Z, int nx, int isgn, int ny, int jsgn);

using namespace std;

class Kvec{
	public:
		double **X,**Y;
		complex<double> **Za;
		int Nx,Ny,ndat;
		double dx[2];
		double Xa[2],Xb[2],Wd[2];
		void set_Xa(double x, double y);
		void set_Wd();
		void set_dx(double x, double y);
		void load(char *fn);
		void load0(char *fn);
		void write(char *fn);
		void write_Za(char fn[128]);
		double freq;
		void clear();
		void init(int n1, int n2);
		void fft(int isgn,int init);
	private:
		double **mem_alloc(int n1,int n2);
	protected:
};
double** Kvec::mem_alloc(int nx, int ny){
	int i;
	double  *p=(double *)malloc(sizeof(double)*nx*ny);
	double **P=(double **)malloc(sizeof(double*)*nx);
	for(i=0;i<nx*ny;i++) p[i]=0.0;
	for(i=0;i<nx;i++) P[i]=p+ny*i;
	return(P);
};
void Kvec::init(int nx,int ny){
	Nx=nx;
	Ny=ny;
	ndat=Nx*Ny;

	int i,j;
	double *p,*q;
	complex<double> *z;
	p=(double *)malloc(sizeof(double)*ndat);
	q=(double *)malloc(sizeof(double)*ndat);
	z=(complex<double> *)malloc(sizeof(complex<double>)*ndat);

	for(i=0;i<ndat;i++){
	       p[i]=0.0;
	       q[i]=0.0;
	       z[i]=complex<double>(0.0,0.0);
	}
	X=(double **)malloc(sizeof(double*)*Nx);
	Y=(double **)malloc(sizeof(double*)*Nx);
	Za=(complex<double> **)malloc(sizeof(complex<double>)*Nx);
	for( i=0;i<Nx;i++){
	      	X[i]=p+Ny*i;
		Y[i]=q+Ny*i;
		Za[i]=z+Ny*i;
	}

};
void Kvec::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void Kvec::set_dx(double x, double y){
	dx[0]=x; dx[1]=y;
};
void Kvec::set_Wd(){
	Xb[0]=Xa[0]+(Nx-1)*dx[0];
	Xb[1]=Xa[1]+(Ny-1)*dx[1];
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
};
void Kvec::load0(char *fn){
	FILE *fp=fopen(fn,"r");
	char cbff[128];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&freq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d,%d\n",&Nx,&Ny);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",dx,dx+1);
	fgets(cbff,128,fp);

	Kvec::set_Wd();
	Kvec::init(Nx,Ny);

	int i,j;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		fscanf(fp,"%le,%le\n",X[i]+j,Y[i]+j);
	}
	}
	fclose(fp);
}
void Kvec::load(char *fn){
	FILE *fp=fopen(fn,"r");
	char cbff[128];
	int nx,ny;

	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&freq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d,%d\n",&nx,&ny);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",dx,dx+1);
	fgets(cbff,128,fp);
	Nx=nx; Ny=ny;
	Kvec::set_Wd();
	double dX[2];
	dX[0]=dx[0]; dX[1]=dx[1];

	double **A=Kvec::mem_alloc(nx,ny);	
	double **B=Kvec::mem_alloc(nx,ny);	
	int p,q;
	p=ceil(log2(nx));
	q=ceil(log2(ny));
	printf("p,q=%d %d\n",p,q);
	printf("nx,ny=%d %d\n",nx,ny);
	Nx=int(pow(2,p));
	Ny=int(pow(2,q));
	printf(" Nx, Ny= %d %d\n",Nx,Ny);

	int i,j;
	for(i=0;i<nx;i++){
	for(j=0;j<ny;j++){
		fscanf(fp,"%le,%le\n",A[i]+j,B[i]+j);
	}
	}
	fclose(fp);

	double xx,yy;
	int l,ix,jy;
	double kx[4],ky[4],phi[4];
	double xi[2],et[2];
	Kvec::init(Nx,Ny);
	dx[0]=(Xb[0]-Xa[0])/(Nx-1);
	dx[1]=(Xb[1]-Xa[1])/(Ny-1);
	printf("dx=%lf %lf\n",dx[0],dx[1]);
	printf("dX=%lf %lf\n",dX[0],dX[1]);
	printf("Wd=%lf %lf\n",Wd[0],Wd[1]);
	double eps=1.e-08;
	//dX[0]*=(1.+eps);
	//dX[1]*=(1.+eps);
	for(i=0; i<Nx; i++){
		xx=fabs(dx[0]*i);
		ix=floor(xx/fabs(dX[0]));
	for(j=0; j<Ny; j++){
		yy=fabs(dx[1]*j);
		jy=floor(yy/fabs(dX[1]));

		if(ix >= nx-1) ix--;
		if(jy >= ny-1) jy--;
		if(ix < 0) ix++;
		if(jy < 0) jy++;
		kx[0]=A[ix][jy];
		kx[1]=A[ix+1][jy];
		kx[2]=A[ix+1][jy+1];
		kx[3]=A[ix][jy+1];

		ky[0]=B[ix][jy];
		ky[1]=B[ix+1][jy];
		ky[2]=B[ix+1][jy+1];
		ky[3]=B[ix][jy+1];

		et[0]=2.*(fabs(xx/dX[0])-ix)-1.0;
		et[1]=2.*(fabs(yy/dX[1])-jy)-1.0;
		phi[0]=0.25*(1.0-et[0])*(1.0-et[1]);
		phi[1]=0.25*(1.0+et[0])*(1.0-et[1]);
		phi[2]=0.25*(1.0+et[0])*(1.0+et[1]);
		phi[3]=0.25*(1.0-et[0])*(1.0+et[1]);
		xi[0]=0.0;
		xi[1]=0.0;
		for(l=0;l<4;l++){
			xi[0]+=(phi[l]*kx[l]);
			xi[1]+=(phi[l]*ky[l]);
		}
		X[i][j]=xi[0];
		Y[i][j]=xi[1];
	}
	}
}
double angle(double x, double y){
	double PI=4.0*atan(1.0);
	double PI2=PI*0.5;
	double v=sqrt(x*x+y*y);
	double alph=asin(y/v);
	if(x<0.0) alph=-PI-alph;
	alph+=PI2;
	return(alph);
}
void Kvec::fft(int isgn,int init){
	int i,j;
	double PI=4.0*atan(1.0);
	double PI2=PI*0.5;
	if(init==1){
		double alph,v;
		for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			/*
			v=X[i][j]*X[i][j]+Y[i][j]*Y[i][j];
			v=sqrt(v);
			alph=asin(Y[i][j]/v);
			if(X[i][j]<0.0) alph=-PI-alph;
			alph+=PI2;
			*/
			alph=angle(X[i][j],Y[i][j]);
		       	Za[i][j]=complex<double>(alph,0.0);
		}
		}
	}else{

		std::mt19937 mt(11);
		std::uniform_real_distribution<double> RndR(0.0,2.*PI);
		double th;
		for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			//Za[i][j]=abs(Za[i][j]*Za[i][j]);
			th=RndR(mt);
			Za[i][j]=complex<double>(cos(th),sin(th))*abs(Za[i][j]);
		}
		}
	};
	FFT2D(Za,Nx,isgn,Ny,isgn);
};
void FFT2D( complex<double> **Z, int nx, int isgn, int ny, int jsgn){
	int i,j;
	int ndat=nx*ny;
	complex<double> *z;
	complex<double> **W; // ny-by-nx complex array
	z=(complex<double> *)malloc(sizeof(complex<double>)*ndat);
	W=(complex<double> **)malloc(sizeof(complex<double>)*ny);
	for(i=0;i<ndat;i++) z[i]=complex<double>(0.0,0.0);
	for(j=0;j<ny;j++) W[j]=z+j*nx;

	for(j=0;j<ny;j++){
		for(i=0;i<nx;i++) W[j][i]=Z[i][j];
		fft(W[j],nx,isgn);
		for(i=0;i<nx;i++) Z[i][j]=W[j][i];
	}

	for(i=0;i<nx;i++) fft(Z[i],ny,jsgn);
};
void Kvec::write(char fn[128]){
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
		fprintf(fp,"%le,%le\n",X[i][j],Y[i][j]);
	}
	}
	fclose(fp);
};
void Kvec::write_Za(char fn[128]){
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
	double ax,ay,alph;
	double PI=4.0*atan(1.0);
	for( i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		ax=Za[i][j].real();
		ay=Za[i][j].imag();
		ax=ax*180./PI;
		ay=ay*180./PI;
		alph=angle(X[i][j],Y[i][j]);
		fprintf(fp,"%le,%le\n",ax,ay);
	}
	}
	fclose(fp);
};


int main(){
	char fname[128]="k117.out";

	Kvec Kx;
	//Kx.load0(fname);
	Kx.load(fname);
	int ifwd,init;
	
	ifwd=1; init=1;
	Kx.fft(ifwd,init);

	ifwd=-1; init=0;
	Kx.fft(ifwd,init);
	char fnout[128]="tmp.dat";
	Kx.write(fnout);
	sprintf(fnout,"tmp2.dat");
	Kx.write_Za(fnout);
	exit(-1);

	double xp[2],xi[2];
	int i0,j0,i,j;
	i0=0;
	j0=20;
	int k=0,kmax=1000;
	double kx[4],ky[4],phi[4],kk;
	double et[2];
	double alph=0.20;
	double ds=Kx.dx[0]*alph;
	int ip,Np=Kx.Nx*5;
	double xp0,yp0;
	double dx[2],Xa[2];

	Xa[0]=Kx.Xa[0]; Xa[1]=Kx.Xa[1];
	dx[0]=Kx.dx[0]; dx[1]=Kx.dx[1];
	xp0=dx[0]*(i0+0.5)+Xa[0];
	yp0=dx[1]*(j0+0.5)+Xa[1];
	int l;
	for(ip=0;ip<Np;ip++){
		xp[0]=xp0;
		xp[1]=yp0;
		printf("%lf %lf\n",xp[0],xp[1]);
	for(k=0;k<kmax;k++){
		i=floor((xp[0]-Xa[0])/dx[0]);
		j=floor((xp[1]-Xa[1])/dx[1]);

		if(i<0) break;
		if(j<0) break;
		if(i>=Kx.Nx-1) break;
		if(j>=Kx.Ny-1) break;
		kx[0]=Kx.X[i][j];
		kx[1]=Kx.X[i+1][j];
		kx[2]=Kx.X[i+1][j+1];
		kx[3]=Kx.X[i][j+1];

		ky[0]=Kx.Y[i][j];
		ky[1]=Kx.Y[i+1][j];
		ky[2]=Kx.Y[i+1][j+1];
		ky[3]=Kx.Y[i][j+1];

		et[0]=2.*((xp[0]-Xa[0])/dx[0]-i)-1.0;
		et[1]=2.*((xp[1]-Xa[1])/dx[1]-j)-1.0;
		phi[0]=0.25*(1.0-et[0])*(1.0-et[1]);
		phi[1]=0.25*(1.0+et[0])*(1.0-et[1]);
		phi[2]=0.25*(1.0+et[0])*(1.0+et[1]);
		phi[3]=0.25*(1.0-et[0])*(1.0+et[1]);
		xi[0]=0.0;
		xi[1]=0.0;
		//phi[0]=1.0;phi[1]=0.0;phi[2]=0.0;phi[3]=0.0;
		for(l=0;l<4;l++){
			xi[0]+=(phi[l]*kx[l]);
			xi[1]+=(phi[l]*ky[l]);
		}
		kk=sqrt(xi[0]*xi[0]+xi[1]*xi[1]);
		xi[0]/=kk;
		xi[1]/=kk;
		xp[0]+=(xi[0]*abs(ds));
		xp[1]+=(xi[1]*abs(ds));
		printf("%lf %lf\n",xp[0],xp[1]);
	}
		xp0+=ds;
		printf("#\n");
	}
	return(0);
};
