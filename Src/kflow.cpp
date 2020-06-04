#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "waves.h"
//#include "kvecs.h"

class Kvec{
	public:
		double **X,**Y;
		int Nx,Ny,ndat;
		double dx[2];
		double Xa[2],Xb[2],Wd[2];
		void set_Xa(double x, double y);
		void set_Wd();
		void set_dx(double x, double y);
		void load(char *fn);
		void load0(char *fn);
		void write(char *fn);
		double freq;
		void clear();
		void init(int n1, int n2);
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
	p=(double *)malloc(sizeof(double)*ndat);
	q=(double *)malloc(sizeof(double)*ndat);

	for(i=0;i<ndat;i++){
	       p[i]=0.0;
	       q[i]=0.0;
	}
	X=(double **)malloc(sizeof(double*)*Nx);
	Y=(double **)malloc(sizeof(double*)*Nx);
	for( i=0;i<Nx;i++) X[i]=p+Ny*i;
	for( i=0;i<Nx;i++) Y[i]=q+Ny*i;

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
	dX[0]*=(1.+eps);
	dX[1]*=(1.+eps);
	for(i=0; i<Nx; i++){
		xx=fabs(dx[0]*i);
		ix=floor(xx/fabs(dX[0]));
	for(j=0; j<Ny; j++){
		yy=fabs(dx[1]*j);
		jy=floor(yy/fabs(dX[1]));

		//if(ix >= nx-1) ix--;
		//if(jy >= ny-1) jy--;
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

int main(){
	char fname[128]="k201.out";

	Kvec Kx;
	//Kx.load0(fname);
	Kx.load(fname);
	char fnout[128]="tmp.dat";
	Kx.write(fnout);
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
