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

void Array2D::out(char *fn){
	FILE *fp=fopen(fn,"w");

	int j;
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
//------------------------------------------------------------
void Array3D::load(char *dir_name){
	FILE *fp;
	char fname[128];

	int i,j,k;

	//sprintf(fname,"%s/scope_%d.csv",dir_name,0);
	//awv.load(fname);

	k=0;
	for(j=0;j<Ny;j++){
	printf("j=%d\n",j);
	bool count=true;
	for(i=0;i<Nx;i++){
		sprintf(fname,"%s/scope_%d.csv",dir_name,k);
		awv.amp=A[i][j];
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

	//int N[3]={Nx,Ny,Nz};
	//printf("Xa=%lf %lf %lf\n",Xa[0],Xa[1],Xa[2]);
	//printf("dx=%lf %lf %lf\n",dx[0],dx[1],dx[2]);
	for(int i=0;i<3;i++){
		Xb[i]=Xa[i]+(Nd[i]-1)*dx[i];
		Wd[i]=Xb[i]-Xa[i];
	}
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
Array3Dcmplx::Array3Dcmplx(int nx, int ny, int nz){
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

#if DB ==1
int main(){

	Grid Gd;	// Measurement Grid 
	char fname[128]="../W20H30_fine/xyl.csv";
	Gd.load(fname); // import grid info.


	int i;
/*
	for(i=5;i+41<=2500;i+=41){
	sprintf(fname,"../W20H30_fine/scope_%d.csv",i);
	Wv1D awv1(fname); // import wave data 1

	sprintf(fname,"../W20H30_fine/scope_%d.csv",i+41);
	Wv1D awv2(fname); // import wave data 2
	Wv1D cor;	// correlation function

	char fn[128];
	cor=corr(awv1,awv2);
	strcpy(fn,"cor.out");
	cor.out_Amp(fn,cor.Np/2); 
	//cor.print_info();
	}
*/


	sprintf(fname,"../W20H30_fine/scope_%d.csv",0);
	Wv1D awv1(fname); // import wave data 1
	Array3D WV(Gd.Nx,Gd.Ny,awv1.Nt);
	WV.print_dim();
	char dir_name[128]="../W20H30_fine";
	WV.load(dir_name);
	WV.awv.print_info();

	//WV.CorrY();
	//WV.CorrX();
	WV.Lp(2);

	return(0);
};
#endif
#if DB ==2
int main(){

	Grid Gd;	// Measurement Grid 
	char fname[128]="../CoreM_short3/x30y20/xyl.csv";
	Gd.load(fname); // import grid info.

	//---------------------------------------
	sprintf(fname,"../CoreM_short3/x30y20/scope_%d.csv",0);
	Wv1D awv1; // import wave data 1
	awv1.load2(fname,true);
	int Nf=awv1.FFT(1);
	Array3D WV(Gd.Nx,Gd.Ny,awv1.Nt);
	Array3Dcmplx WVf(Gd.Nx,Gd.Ny,Nf);
	Array3D Phi(Gd.Nx,Gd.Ny,Nf/2);
	Array2D Kx(Gd.Nx,Gd.Ny);
	Array2D Ky(Gd.Nx,Gd.Ny);

	WV.set_Xa(Gd.Xcod[0],Gd.Ycod[0],awv1.t1);
	WV.set_dx(Gd.dx,Gd.dy,awv1.dt);
	WV.set_Wd();

	WV.print_dim();
	WVf.print_dim();
	char dir_name[128]="../CoreM_short3/x30y20";
	WV.load(dir_name);

	int i,j,k;
	for(i=0;i<Gd.Nx;i++){
		printf("i=%d/%d\n",i,Gd.Nx);
	for(j=0;j<Gd.Ny;j++){
		awv1.amp=WV.A[i][j];
		awv1.Amp=WVf.Z[i][j];
		awv1.FFT(1);
	}
	}
	char fout[128]="wvf.out";
	WVf.write_zslice(fout,int(1.0/awv1.df));
	for(i=0;i<Gd.Nx;i++){
	for(j=0;j<Gd.Ny;j++){
	for(k=0;k<Nf/2;k++){
		Phi.A[i][j][k]=arg(WVf.Z[i][j][k]);
	}
	}
	}
	k=int(1.2/awv1.df);
	double dpx,dpy,wgt;
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);
	wgt=2.0*Gd.dx;
	for(i=0;i<Gd.Nx-1;i++){
	for(j=0;j<Gd.Ny;j++){
		dpx=Phi.A[i+1][j][k]-Phi.A[i][j][k];
		if(dpx> PI) dpx=-(PI2-dpx);
		if(dpx<-PI) dpx=PI2+dpx;
		dpx/=wgt;
		Kx.A[i][j]+=dpx;
		Kx.A[i+1][j]+=dpx;
	}
	}
	for(j=0;j<Gd.Ny;j++){
		Kx.A[0][j]*=2.0;
		Kx.A[Gd.Nx-1][j]*=2.0;
	}

	wgt=2.0*Gd.dy;
	for(i=0;i<Gd.Nx;i++){
	for(j=0;j<Gd.Ny-1;j++){
		dpy=Phi.A[i][j+1][k]-Phi.A[i][j][k];
		if(dpy> PI) dpy=-(PI2-dpy);
		if(dpy<-PI) dpy=PI2+dpy;
		dpy/=wgt;
		Ky.A[i][j]+=dpy;
		Ky.A[i][j+1]+=dpy;
	}
	}
	for(i=0;i<Gd.Nx;i++){
		Ky.A[i][0]*=2.0;
		Ky.A[i][Gd.Ny-1]*=2.0;
	}
	char fnkx[128]="kx.out";
	char fnky[128]="ky.out";
	Kx.out(fnkx);
	Ky.out(fnky);


	exit(-1);
	//WV.Butterworth(0.0,3.0);

	//---------------------------------------
	Array2D Bwv,Tg;
	Bwv=WV.proj();
	puts("done");
	char fn[128]="bwv.out";
	Bwv.out(fn);

	Tg=WV.gdelay(2.9);
	sprintf(fn,"tg.out");
	Tg.out(fn);
	
	return(0);
};
#endif
