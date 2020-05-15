#define DB 3
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
//#include "fft.h"
//
#include "waves.h"

using namespace std;

//------------------------------------------------------------
int Grid::load(char *fname){
	FILE *fp;
	char cbff[128];

	fp=fopen(fname,"r");
	int ndat=0;
	while(fgets(cbff,128,fp) != NULL) ndat++;
	fclose(fp);

	printf("ndat=%d\n",ndat);

	Xcod=(double *)malloc(sizeof(double)*ndat);
	Ycod=(double *)malloc(sizeof(double)*ndat);
	vals=(int *)malloc(sizeof(int)*ndat);
	fp=fopen(fname,"r");
	int i;
	for(i=0;i<ndat;i++){
		fscanf(fp,"%lf,%lf,%d\n",Xcod+i,Ycod+i,vals+i);
	}
	
	double y0=Ycod[0];
	Nx=1;
	while(Ycod[Nx]==y0){
		Nx++;
	}
	Ny=ndat/Nx;
	dx=Xcod[1]-Xcod[0];
	dy=Ycod[Nx]-Ycod[0];
	printf("Nx=%d, Ny=%d\n",Nx,Ny);
	printf("dx=%lf, dy=%lf\n",dx,dy);
	printf("Ycod[%d]=%lf, Ycod[0]=%lf\n",Nx,Ycod[Nx],Ycod[0]);
	
	fclose(fp);
	return(ndat);
};
//------------------------------------------------------------
Array2D::Array2D(int nx,int ny){
	Nx=nx;
	Ny=ny;
	ndat=Nx*Ny;

	int i,j;
	A2=(double *)malloc(sizeof(double)*ndat);
	for(i=0;i<ndat;i++) A2[i]=0.0;

	A=(double **)malloc(sizeof(double*)*Nx);
	for( i=0;i<Nx;i++) A[i]=A2+Ny*i;

	for(i=0;i<2;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}

};
/*
Array2D::Array2D(){
	Nx=1;
	Ny=1;
	ndat=Nx*Ny;
	A2=(double *)malloc(sizeof(double)*ndat);
	A=(double **)malloc(sizeof(double*)*Nx);
	for(int i=0;i<2;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}
};
*/
Array2D::Array2D(){};
void Array2D::mem_alloc(){
	int i;
	ndat=Nx*Ny;
	A2=(double *)malloc(sizeof(double)*ndat);
	for(i=0;i<ndat;i++) A2[i]=0.0;
	A=(double **)malloc(sizeof(double*)*Nx);
	for(i=0;i<Nx;i++) A[i]=A2+Ny*i;
};
/*
Array2D::~Array2D(){
	free(A);
	free(A2);
};
*/
void Array2D::set_Xa(double x, double y){
	Xa[0]=x; Xa[1]=y;
};
void Array2D::set_dx(double x, double y){
	dx[0]=x; dx[1]=y;
};
void Array2D::set_Wd(){

	Xb[0]=Xa[0]+(Nx-1)*dx[0];
	Xb[1]=Xa[1]+(Ny-1)*dx[1];
	Wd[0]=Xb[0]-Xa[0];
	Wd[1]=Xb[1]-Xa[1];
};

void Array2D::load(char *fn){
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

	Array2D::set_Wd();
	Array2D::mem_alloc();

	int i,j;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		fscanf(fp,"%le\n",A[i]+j);
	}
	}

	fclose(fp);
}
void Array2D::out(char *fn){
	FILE *fp=fopen(fn,"w");

	int j;
	fprintf(fp,"# frequency [MHz]\n");
	fprintf(fp,"%lf\n",freq);
	fprintf(fp,"# Nx, Ny\n");
	fprintf(fp,"%d,%d\n",Nx,Ny);
	fprintf(fp,"# Xa[0:1]\n");
	fprintf(fp,"%lf,%lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# dx[0:1]\n");
	fprintf(fp,"%lf,%lf\n",dx[0],dx[1]);
	fprintf(fp,"# amp (for x{ for y})\n");
	for( int i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fprintf(fp,"%le\n",A[i][j]);
	}
	}
	fclose(fp);
};
void Array2D::clear(){
	int i,j;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		A[i][j]=0.0;
	}
	}
};
//------------------------------------------------------------
void Array3D::load(char *dir_name){
	FILE *fp;
	char fname[128];

	int i,j,k;

	//sprintf(fname,"%s/scope_%d.csv",dir_name,0);
	//awv.load(fname);

	k=0;
	bool count=true;
	int nn;
	for(j=0;j<Ny;j++){
		printf("j=%d\n",j);
	for(i=0;i<Nx;i++){
		sprintf(fname,"%s/scope_%d.csv",dir_name,k);
		awv.amp=A[i][j];
		awv.mllc=true;
		awv.load2(fname,count);
		k++;
		count=false;
	}
	}

	//for(j=0;j<awv.Nt;j++) printf("%le\n",A[0][0][j]);
};
Array3D::Array3D(int nx, int ny, int nz){
	Nx=nx;
	Ny=ny;
	Nz=nz;

	Nd[0]=Nx;
	Nd[1]=Ny;
	Nd[2]=Nz;
	
	ndat=Nx*Ny*Nz;
	int i,j,k;  

	A3=(double *)malloc(sizeof(double)*ndat);
	for(k=0;k<ndat;k++) A3[k]=0.0;

	A2=(double **)malloc(sizeof(double *)*nx*ny);	
	for(j=0;j<nx*ny; j++) A2[j]=A3+j*nz;

	A=(double ***)malloc(sizeof(double **)*nx);
	for(i=0;i<nx;i++) A[i]=A2+i*ny;

	for(i=0;i<3;i++){
		dx[i]=1.0;
		Xa[i]=0.0;
	}
}
void Array3D::out(char *fname){
	FILE *fp=fopen(fname,"w");
	fprintf(fp,"# Xa[0:2]\n");
	fprintf(fp,"%lf, %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# dx[0:2]\n");
	fprintf(fp,"%lf, %lf\n",dx[0],dx[1]);
	fprintf(fp,"# Nd[0:2]\n");
	fprintf(fp,"%d, %d\n",Nd[0],Nd[1]);
	fprintf(fp,"# t1, dt, Nt\n");
	fprintf(fp,"%lf, %lf, %d\n",Xa[2],dx[2],Nd[2]);
	fprintf(fp,"# amp\n");
	int i,j,k;
	for(i=0;i<Nd[0];i++){
	for(j=0;j<Nd[1];j++){
	for(k=0;k<Nd[2];k++){
		fprintf(fp,"%le\n",A[i][j][k]);
	}
	}
	}
	fclose(fp);
}
Array3D::~Array3D(){
	free(A3);
	free(A2);
	free(A);
};
void Array3D::set_Xa(double x, double y, double z){
	Xa[0]=x; Xa[1]=y; Xa[2]=z;
};
void Array3D::set_dx(double x, double y, double z){
	dx[0]=x; dx[1]=y; dx[2]=z;
};
void Array3D::set_Wd(){
	for(int i=0;i<3;i++){
		Xb[i]=Xa[i]+(Nd[i]-1)*dx[i];
		Wd[i]=Xb[i]-Xa[i];
	}
};
int Array3D::get_index(double val, int axis){
	int indx=round((val-Xa[axis])/dx[axis]);
	if(indx <0) indx=-1;
	if(indx >=Nd[axis]) indx=-1;
	return(-1);
};
void Array3D::print_dim(){
	printf("Array size=(%d, %d, %d)\n",Nx,Ny,Nz);
};
void Array3D::set_val(double v){
	for(int i=0;i<ndat;i++) A3[i]=v;
};
void Array3D::Lp(int p){
	int i,j;
	Array2D Lp(Nx,Ny);

	Wv1D wv;
	wv=awv;
	wv.print_info();

	double t1=wv.t1;
	double t2=wv.t2;
	char fn[128];
	if(p==2){
		sprintf(fn,"L2.out");
		for(i=0;i<Nx;i++){
			printf("i=%d/%d\n",i+1,Nx);
		for(j=0;j<Ny;j++){
			wv.amp=A[i][j];
			Lp.A[i][j]=wv.L2(t1,t2);
		}
		}
	}else if(p==0){
		sprintf(fn,"Linf.out");
		for(i=0;i<Nx;i++){
			printf("i=%d/%d\n",i+1,Nx);
		for(j=0;j<Ny;j++){
			wv.amp=A[i][j];
			Lp.A[i][j]=wv.max(t1,t2);
		}
		}
	};
	Lp.out(fn);
};
Array2D Array3D::proj(){
	Array2D Bdat(Ny,Nz);
	int i,j,k;
	//printf("Nx,Ny,Nz=%d %d %d\n",Nx,Ny,Nz);
	for(i=0; i<Ny; i++){
		//printf("i=%d\n",i);
	for(j=0; j<Nz; j++){
		for(k=0; k<Nx; k++) Bdat.A[i][j]+=A[k][i][j]; 
		Bdat.A[i][j]/=Nx;
	}
	}
	Bdat.set_Xa(Xa[1],Xa[2]);
	Bdat.set_dx(dx[1],dx[2]);
	Bdat.set_Wd();
	return(Bdat);
};
void Array3D::CorrY(){

	int i,j;
	double tmax,Amax;
	Array2D Tm(Nx,Ny-1);
	Array2D Amp(Nx,Ny-1);
	char fn[128]="debug.dat";

	Wv1D wv1,wv2,cor;
	wv1=awv;
	wv2=awv;
	wv1.print_info();
	wv2.print_info();
	int tmp;
	for(i=0;i<Nx;i++){
		printf("i=%d/%d\n",i+1,Nx);
	for(j=0;j<Ny-1;j++){
		wv1.amp=A[i][j];
		wv2.amp=A[i][j+1];
		cor=corr(wv1,wv2,&tmax,&Amax);
		Tm.A[i][j]=tmax;
		Amp.A[i][j]=Amax;
	}
	}

	char fname[128];
	strcpy(fname,"tmax.out");
	Tm.out(fname);
	strcpy(fname,"amax.out");
	Amp.out(fname);

};
void Array3D::CorrX(){

	int i,j;
	double tmax,Amax;
	Array2D Tm(Nx-2,Ny);
	Array2D Amp(Nx-2,Ny);

	Wv1D wv1,wv2,cor;
	wv1=awv;
	wv2=awv;
	wv1.print_info();
	wv2.print_info();
	int tmp;
	for(i=0;i<Nx-2;i++){
		printf("i=%d/%d\n",i+1,Nx-1);
	for(j=0;j<Ny;j++){
		wv1.amp=A[i][j];
		wv2.amp=A[i+2][j];
		cor=corr(wv1,wv2,&tmax,&Amax);
		Tm.A[i][j]=tmax;
		Amp.A[i][j]=Amax;
	}
	}

	char fname[128];
	strcpy(fname,"tmax2.out");
	Tm.out(fname);
	strcpy(fname,"amax2.out");
	Amp.out(fname);
};
void Array3D::Butterworth(double cx, double cy){
	int i,j;
	double x,y;
	double tb;
	double Tw_6dB=5.0,dlt=12.5;
	for(i=0;i<Nx;i++){
		x=Xa[0]+dx[0]*i;
	for(j=0;j<Ny;j++){
		y=Xa[1]+dx[1]*j;
		tb=y/cy+dlt+Tw_6dB*0.5;
		awv.amp=A[i][j];
		awv.Butterworth(tb,Tw_6dB);
	}
	}

};
Array2D Array3D::gdelay( double cy){
	Array2D Tg(Nx,Ny); 
	Tg.set_Xa(Xa[0],Xa[1]);
	Tg.set_dx(dx[0],dx[1]);
	Tg.set_Wd();
	int i,j;
	double x,y,tb,tgb;
	double Tw_6dB=3.0,dlt=12.5;
	for(i=0;i<Nx;i++){
		printf("i=%d\n",i);
		x=Xa[0]+dx[0]*i;
	for(j=0;j<Ny;j++){
		y=Xa[1]+dx[1]*j;
		tb=y/cy+dlt+Tw_6dB*0.5;
		awv.amp=A[i][j];
		awv.Butterworth(tb,Tw_6dB);
		Tg.A[i][j]=awv.gdelay();
		awv.fft_stat=0;
	}
	}
	puts("done");
	return(Tg);
};
//-------------------------------------------------------
Array3Dcmplx::Array3Dcmplx(int nx, int ny, int nz){
	/*
	Nx=nx;
	Ny=ny;
	Nz=nz;

	Nd[0]=Nx;
	Nd[1]=Ny;
	Nd[2]=Nz;
	
	ndat=Nx*Ny*Nz;
	int i,j,k;  

	Z3=(complex<double> *)malloc(sizeof(complex<double>)*ndat);
	for(k=0;k<ndat;k++) Z3[k]=complex<double>(0.0,0.0);

	Z2=(complex<double> **)malloc(sizeof(complex<double> *)*nx*ny);	
	for(j=0;j<nx*ny; j++) Z2[j]=Z3+j*nz;

	Z=(complex<double> ***)malloc(sizeof(complex<double> **)*nx);
	for(i=0;i<nx;i++) Z[i]=Z2+i*ny;
	*/

	Array3Dcmplx::mem_alloc(nx,ny,nz);
};
Array3Dcmplx::Array3Dcmplx(){};
void Array3Dcmplx::mem_alloc(int nx, int ny, int nz){
	Nx=nx;
	Ny=ny;
	Nz=nz;

	Nd[0]=Nx;
	Nd[1]=Ny;
	Nd[2]=Nz;
	
	ndat=Nx*Ny*Nz;
	int i,j,k;  

	Z3=(complex<double> *)malloc(sizeof(complex<double>)*ndat);
	for(k=0;k<ndat;k++) Z3[k]=complex<double>(0.0,0.0);

	Z2=(complex<double> **)malloc(sizeof(complex<double> *)*nx*ny);	
	for(j=0;j<nx*ny; j++) Z2[j]=Z3+j*nz;

	Z=(complex<double> ***)malloc(sizeof(complex<double> **)*nx);
	for(i=0;i<nx;i++) Z[i]=Z2+i*ny;
};
int Array3Dcmplx::get_index(double val, int axis){
	int indx=round((val-Xa[axis])/dx[axis]);
	if(indx <0) indx=-1;
	if(indx >=Nd[axis]) indx=-1;
	return(indx);
};

void Array3Dcmplx::out(char *fname){
	FILE *fp=fopen(fname,"w");
	fprintf(fp,"# Xa[0:2]\n");
	fprintf(fp,"%lf, %lf\n",Xa[0],Xa[1]);
	fprintf(fp,"# dx[0:2]\n");
	fprintf(fp,"%lf, %lf\n",dx[0],dx[1]);
	fprintf(fp,"# Nd[0:2]\n");
	fprintf(fp,"%d, %d\n",Nd[0],Nd[1]);
	fprintf(fp,"# f1, df, Nf\n");
	fprintf(fp,"%lf, %lf, %d\n",Xa[2],dx[2],Nd[2]);
	fprintf(fp,"# amp(Re, Im)\n");
	int i,j,k;
	for(i=0;i<Nd[0];i++){
	for(j=0;j<Nd[1];j++){
	for(k=0;k<Nd[2]/2;k++){
		fprintf(fp,"%le, %le\n",real(Z[i][j][k]),imag(Z[i][j][k]));
	}
	}
	}
	fclose(fp);
};
void Array3Dcmplx::load(char *fname){
	char cbff[128];
	FILE *fp=fopen(fname,"r");
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf\n",dx,dx+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%d, %d\n",Nd,Nd+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf, %lf, %d\n",Xa+2,dx+2,Nd+2);
	fgets(cbff,128,fp);

	Array3Dcmplx::mem_alloc(Nd[0],Nd[1],Nd[2]);
	int i,j,k;
	double vr,vi;
	for(i=0;i<Nd[0];i++){
	for(j=0;j<Nd[1];j++){
	for(k=0;k<Nd[2]/2;k++){
		fscanf(fp,"%le, %le\n",&vr,&vi);
		Z[i][j][k]=complex<double>(vr,vi);
	}
	}
	}
	fclose(fp);
};
int Array3Dcmplx::write_zslice(char *fname, int k){

	FILE *fp=fopen(fname,"w");
	if(k>=Nz) return(0);

	int i,j;
	fprintf(fp,"%d,%d\n",Nx,Ny);
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		fprintf(fp,"%lf, %lf\n",Z[i][j][k].real(),Z[i][j][k].imag());
	}
	}
	fclose(fp);
	return(Nx*Ny);

};
Array3Dcmplx::~Array3Dcmplx(){
	free(Z3);
	free(Z2);
	free(Z);
};
void Array3Dcmplx::print_dim(){
	printf("Array size=(%d, %d, %d)\n",Nx,Ny,Nz);
};
void Array3Dcmplx::set_Xa(double x, double y, double z){
	Xa[0]=x; Xa[1]=y; Xa[2]=z;
};
void Array3Dcmplx::set_dx(double x, double y, double z){
	dx[0]=x; dx[1]=y; dx[2]=z;
};
void Array3Dcmplx::set_Wd(){
	for(int i=0;i<3;i++){
		Xb[i]=Xa[i]+(Nd[i]-1)*dx[i];
		Wd[i]=Xb[i]-Xa[i];
	}
};
void Array3Dcmplx::set_val(double vr, double vi){
	for(int i=0;i<ndat;i++) Z3[i]=complex<double>(vr,vi);
};

