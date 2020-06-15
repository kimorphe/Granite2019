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

int moments(double *dat, int ndat, double *m0, double *m2, double *min, double *max){

	int i,isum=0;
	double M0,M2,Min,Max;
	M0=0.0;
	M2=0.0;
	Min=-1.0;
	Max=-1.0;
	for(i=0;i<ndat;i++){
		if(dat[i]<0.0) continue;
		M0+=dat[i];
		M2+=(dat[i]*dat[i]);
		if(isum==0){
			Min=dat[i];
			Max=dat[i];
		}
		isum++;
		if(Min > dat[i]) Min=dat[i];
		if(Max < dat[i]) Max=dat[i];
	}
	if(isum>0){
		M2=M0*M0-M2;
		M0/=isum;
		M2/=isum;
	}else{
		M0=-1.0;
		M2=-1.0;
	};

	*m0=M0;
	*m2=M2;
	*min=Min;
	*max=Max;
	return(isum);
};
void xproj_mean(double **A,  double *ay, int *ncnt, double thr, int nx, int ny){

	int i,j,k;
	int count;
	for(j=0;j<ny;j++){
		count=0;
		ay[j]=0.0;
	for(i=0;i<nx;i++){
		if(A[i][j]<thr) continue;
		ay[j]+=A[i][j];
		count++;
	}
		if(count==0){
			ay[j]=thr-1.0;
		}else{
			ay[j]/=count;
		}
		ncnt[j]=count;
	}
};
void xproj_min(double **A,  double *ay, int *ncnt, double thr, int nx, int ny){

	int i,j,k;
	int count;
	double amin;
	for(j=0;j<ny;j++){
		count=0;
		ay[j]=0.0;
	for(i=0;i<nx;i++){
		if(A[i][j]<thr) continue;
		if(count==0) amin=A[i][j];

		if(amin> A[i][j]) amin=A[i][j];
		count++;
	}
		if(count==0){
			ay[j]=thr-1.0;
		}else{
			ay[j]=amin;
		}
		ncnt[j]=count;
	}
}

double** FSLICE::mem_alloc(int n, int m){
	int i;
	double **A,*A2;
	A2=(double *)malloc(sizeof(double)*n*m);
	A=(double **)malloc(sizeof(double*)*n);
	for(i=0;i<n*m;i++) A2[i]=0.0;
	for(i=0;i<n;i++)  A[i]=A2+m*i;
	return(A);
};
double*** FSLICE::mem_alloc3d(int n, int m, int l){
	int i;
	double ***A,**A2, *A3;
	A3=(double  *)malloc(sizeof(double)*n*m*l);
	A2=(double **)malloc(sizeof(double*)*n*m);
	A=(double ***)malloc(sizeof(double**)*n);
	for(i=0;i<n*m*l;i++)	A3[i]=0.0;
	for(i=0;i<n*m;i++)	A2[i]=A3+l*i;
	for(i=0;i<n;i++)	 A[i]=A2+m*i;
	return(A);
};
int** FSLICE::mem_ialloc(int n, int m){
	int i;
	int **A,*A2;
	A2=(int *)malloc(sizeof(int)*n*m);
	A=(int **)malloc(sizeof(int*)*n);
	for(i=0;i<n*m;i++) A2[i]=0;
	for(i=0;i<n;i++)  A[i]=A2+m*i;
	return(A);
};

void FSLICE::init(int n, int m){
	Nx=n; Ny=m; ndat=n*m;
	Nd[0]=n; Nd[1]=m;
	malloc_arrays();
};

void FSLICE::malloc_arrays(){
	Kx=mem_alloc(Nx,Ny);
	Ky=mem_alloc(Nx,Ny);
	Kxx=mem_alloc(Nx,Ny);
	Kyy=mem_alloc(Nx,Ny);
	Kxy=mem_alloc(Nx,Ny);
	Phi=mem_alloc(Nx,Ny);
	psi=mem_alloc(Nx,Ny);
	Amp=mem_alloc(Nx,Ny);
	Pave=mem_alloc(Nx,Ny);	// phase (mean) 
	Pvar=mem_alloc(Nx,Ny);	// phase (var)
	Pmin=mem_alloc(Nx,Ny);	// phase (min.)
	Pmax=mem_alloc(Nx,Ny);	// phase (max.)
	Ppnt=mem_alloc(Nx,Ny);	// phase (point source)
	is_alloc=1;
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
double darg(double th1, double th2){
	complex<double> z1=complex<double>(cos(th1),sin(th1));
	complex<double> z2=complex<double>(cos(th2),sin(th2));
	return(arg(z2/z1) );
};
void FSLICE::Integrate2(){ // distance function generation
	Heap hp;
	int i,j,k,kmin;
	int istart;
	double phi0,phi,dphi,phi_tmp;
	/*
	phi0=0.0;
	phi=0.0;
	imin=0;
	for(i=1;i<Nx;i++){
		dphi=darg(Phi[i-1][0],Phi[i][0]);
		phi+=dphi;
		if(phi < phi0){
			imin=i;
			phi0=phi;
		}
	}
	*/

	//imin=50;

	int ix0,jy0, ix,jy;
	int iofst[4]={0,1,0,-1};
	int jofst[4]={-1,0,1,0};
	int **Ix=FSLICE::mem_ialloc(Nx,Ny);
	int nsize=Nx*Ny*0.2;
	int itr,itr_max=Nx*Ny;
	hp.init(nsize);

	double ***Psi=FSLICE::mem_alloc3d(Nx,Ny,Nx);
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	for(k=0;k<Nx;k++){
		Psi[i][j][k]=-1.0;
	}
	}
	}

	for(istart=0;istart<Nx;istart++){
		for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			Ix[i][j]=0;
			psi[i][j]=-1.0;
		}
		}

		ix0=istart; jy0=0;
		kmin=ix0*Ny+jy0;
		hp.ndat=0;
		hp.add(0.0,kmin);
		itr=0;
		while(hp.ndat>0){
			// Pop smallest phase cell
			phi0=hp.pop(&kmin);
			ix0=kmin/Ny;
			jy0=kmin%Ny;
			psi[ix0][jy0]=phi0;
			Ix[ix0][jy0]=1;
			// Create neighbor list
			//for(i=0;i<4;i++){
			for(i=1;i<4;i++){
				ix=ix0+iofst[i];
				if(ix<0) continue;
				if(ix>=Nx) continue;
				jy=jy0+jofst[i];
				if(jy<0) continue;
				if(jy>=Ny) continue;

				if(Ix[ix][jy]==1) continue;
				k=ix*Ny+jy;
				dphi=darg(Phi[ix0][jy0],Phi[ix][jy]);	
				if(dphi >= 0.0){
					phi_tmp=psi[ix0][jy0]+dphi;
					if(Ix[ix][jy]==0){
			       			hp.add(phi_tmp,k);
						Ix[ix][jy]=-1;
					}else if(phi_tmp < psi[ix][jy]){
					       	hp.del(k);
			       			hp.add(phi_tmp,k);
					}
				}
			}
			//printf("phi=%lf, (ix,jy)=(%d,%d) k=%d| ",phi0,ix0,jy0,ix0*Ny+jy0);
			//printf("ndat=%d\n",hp.ndat);
			itr++;
			if(itr>itr_max){
				puts("Too many iteration in Integration2");
				puts("--> abort process.");
				exit(-1);
			}
		}
//		printf("itr=%d\n",itr);
//		char tmp[128]="phix.out";
//		FSLICE::export_phix(tmp);

		for(i=0;i<Nx;i++){
		for(j=0;j<Ny;j++){
			Psi[i][j][istart]=psi[i][j];
		}
		}

	}
	double pmin,isum;
	double m0,m2,min,max;
	int isrc=Nx/2;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		psi[i][j]=-1;
		isum=0;
		pmin=-1;
		Ppnt[i][j]=Psi[i][j][isrc];
		for(k=0;k<Nx;k++){
			if(Psi[i][j][k]<0.0) continue;
			if(isum==0) pmin=Psi[i][j][k];
			isum++;
			if(Psi[i][j][k]<pmin)  pmin=Psi[i][j][k];
		};
		moments(Psi[i][j],Nx,&m0,&m2,&min,&max);
		Pave[i][j]=m0;
		Pvar[i][j]=m2;
		Pmin[i][j]=min;
		Pmax[i][j]=max;
	}
	}
	char tmp[128]="phix.out";
	FSLICE::export_phix(tmp);

	double pb;
	for(j=0;j<Ny;j++){
		pb=0.0;
		isum=0;
	for(i=0;i<Nx;i++){
		if(psi[i][j]<0.0) continue;
		isum++;
		pb+=psi[i][j];
	}
		if(isum>0) pb/=isum;
//		printf("%lf\n",pb);
	}

	free(Psi);
};
void FSLICE::Integrate(){
	int i,j;
	double th1,th2,dpx,dpy;
	complex<double> z1,z2;
	int *A2=(int *)malloc(sizeof(int)*Nx*Ny);
	int **A=(int **)malloc(sizeof(int*)*Nx);
	for(i=0;i<Nx*Ny;i++) A2[i]=0;
	for(i=0;i<Nx;i++) A[i]=A2+Ny*i;


	A[0][0]=1;
	for(i=0;i<Nx-1;i++){
		th1=Phi[i][0];
		th2=Phi[i+1][0];
		z1=complex<double>(cos(th1),sin(th1));
		z2=complex<double>(cos(th2),sin(th2));
		dpx=arg(z2/z1);
		psi[i+1][0]=psi[i][0]+dpx;
		A[i+1][0]=1;
	}

	int nflip;
	for(j=1;j<Ny;j++){

		for(i=0;i<Nx;i++){
			th2=Phi[i][j];
			th1=Phi[i][j-1];
			z1=complex<double>(cos(th1),sin(th1));
			z2=complex<double>(cos(th2),sin(th2));
			dpy=arg(z2/z1);
			A[i][j-1]=abs(A[i][j-1]);
			if(dpy>=0.0){
				if(A[i][j-1]>0){
					psi[i][j]=psi[i][j-1]+dpy;
					A[i][j]=1;
				}
			}
		}

		nflip=1;
		while(nflip>0){
			nflip=0;
			for(i=0;i<Nx-1;i++){
				//A[i][j]=abs(A[i][j]);		
				//A[i+1][j]=abs(A[i+1][j]);		
				if(A[i][j]==0) continue;
				if(A[i+1][j]==1) continue;
				th1=Phi[i][j];
				th2=Phi[i+1][j];
				z1=complex<double>(cos(th1),sin(th1));
				z2=complex<double>(cos(th2),sin(th2));
				dpx=arg(z2/z1);
				if(dpx >=0.0){
				       A[i+1][j]=-1;
				       psi[i+1][j]=psi[i][j]+dpx;
				       nflip++;
				}
			}
			double psid;
			for(i=Nx-1;i>0;i--){
				if(A[i][j]<=0) continue;
				if(A[i-1][j]==1) continue;
				th1=Phi[i][j];
				th2=Phi[i-1][j];
				z1=complex<double>(cos(th1),sin(th1));
				z2=complex<double>(cos(th2),sin(th2));
				dpx=arg(z2/z1);
				if(dpx >= 0.0){
					psid=psi[i][j]+dpx;
					if(A[i-1][j]==0){
						psi[i-1][j]=psid;
						A[i-1][j]=1;
						nflip++;
					}
					if(A[i-1][j]==-1){
						if(psi[i-1][j]>psid){
						      psi[i-1][j]=psid;
						      A[i-1][j]=1;
						      nflip++;
						}
					}
				}
			}
			for(i=0;i<Nx;i++) A[i][j]=abs(A[i][j]);
		}
	}

	char tmp[128]="phix.out";
	FSLICE::export_phix(tmp);
}
void FSLICE::export_phix(char fn[128]){
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
		fprintf(fp,"%le\n",psi[i][j]);
	}
	}
	fclose(fp);
};
void FSLICE::Grad(){
		
	double PI=4.0*atan(1.0);
	double PI2=8.0*atan(1.0);
	double wgt=2.0*dx[0]*PI2;
	double dpx,dpy;
	int i,j;
	for(i=0;i<Nx;i++){ 
	for(j=0;j<Ny;j++){
		Kx[i][j]=0.0; 
		Ky[i][j]=0.0; 
		Kxx[i][j]=0.0;
		Kyy[i][j]=0.0;
		Kxy[i][j]=0.0;
	}
	}

	double dx2=dx[0]*dx[0];
	double dy2=dx[1]*dx[1];

	complex<double> z1,z2;
	double th1,th2;
	for(i=0;i<Nx-1;i++){ 
	for(j=0;j<Ny;j++){
		th1=Phi[i][j];
		th2=Phi[i+1][j];
		z1=complex<double>(cos(th1),sin(th1));
		z2=complex<double>(cos(th2),sin(th2));
		dpx=arg(z2/z1);

		//dpx=Phi[i+1][j]-Phi[i][j];
		//if(dpx> PI) dpx=-(PI2-dpx);
		//if(dpx<-PI) dpx=PI2+dpx;
		Kxx[i][j]+=(dpx/dx2);
		Kxx[i+1][j]-=(dpx/dx2);

		dpx/=wgt;
		Kx[i][j]+=dpx;
		Kx[i+1][j]+=dpx;

		//Kyx[i][j]+=dpx
	}
	}

	for(j=0;j<Ny;j++){
		Kx[0][j]*=2.0;
		Kx[Nx-1][j]*=2.0;

		Kxx[0][j]=0.0;
		Kxx[Nx-1][j]=0.0;
	}

	wgt=2.0*dx[1]*PI2;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny-1;j++){
		th2=Phi[i][j+1];
		th1=Phi[i][j];
		z1=complex<double>(cos(th1),sin(th1));
		z2=complex<double>(cos(th2),sin(th2));
		dpy=arg(z2/z1);

		//dpy=Phi[i][j+1]-Phi[i][j];
		//if(dpy> PI) dpy=-(PI2-dpy);
		//if(dpy<-PI) dpy=PI2+dpy;
		Kyy[i][j]+=(dpy/dy2);
		Kyy[i][j+1]-=(dpy/dy2);

		dpy/=wgt;
		Ky[i][j]+=dpy;
		Ky[i][j+1]+=dpy;
	}
	}
	for(i=0;i<Nx;i++){
		Ky[i][0]*=2.0;
		Ky[i][Ny-1]*=2.0;

		Kyy[i][0]=0.0;
		Kyy[i][Ny-1]=0.0;
	}

	for(i=1;i<Nx-1;i++){
	for(j=1;j<Ny-1;j++){
		Kxy[i][j]=(Kx[i][j+1]-Kx[i][j-1])/(2.*dx[1]);
		Kxy[i][j]=Kxy[i][j]+(Ky[i+1][j]-Ky[i-1][j])/(2.*dx[0]);
		Kxy[i][j]*=0.5;
	}
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
void FSLICE::load_Grad(char fn[128]){
	FILE *fp=fopen(fn,"r");
	int i,j;
	char cbff[128];
	fgets(cbff,128,fp);
	fscanf(fp,"%lf\n",&freq);
	printf("#freq=%lf\n",freq);
	fgets(cbff,128,fp);
	fscanf(fp,"%d,%d\n",&Nx,&Ny);
	printf("#Nx,Ny=%d,%d\n",Nx,Ny);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",Xa,Xa+1);
	fgets(cbff,128,fp);
	fscanf(fp,"%lf,%lf\n",dx,dx+1);
	fprintf(fp,"# kx, ky (for x{ for y})\n");
	fgets(cbff,128,fp);
	double kx,ky;
	if(is_alloc==0){
		Kx=mem_alloc(Nx,Ny);
		Ky=mem_alloc(Nx,Ny);
	};
	for( i=0;i<Nx;i++){
	for( j=0;j<Ny;j++){
		fscanf(fp,"%le,%le\n",&kx, &ky);
		Kx[i][j]=kx;
		Ky[i][j]=ky;
	}
	}
	fclose(fp);
};
void FSLICE::export_Hess(char fn[128]){
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
		fprintf(fp,"%le,%le,%le\n",Kxx[i][j],Kyy[i][j],Kxy[i][j]);
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
	k_mean=0.0; 
	k_sig=0.0; 
	a_mean=0.0; 
	a_sig=0.0;
	asum=0.0;
	for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
		xi0=Kx[i][j];
		xi1=Ky[i][j];
		xi=sqrt(xi0*xi0+xi1*xi1);
		alph=asin(xi1/xi);

		if(xi1>0.0) continue;

		if(xi0<0.0) alph=-PI-alph;
		alph=alph/PI*180.0;

		amp=Amp[i][j];
		amp*=amp;
		asum+=amp;

		//ibin=(xi-kmin)/dk;
		ibin=(abs(xi1)-kmin)/dk;
		//if(ibin>=0 && ibin <nbin) Prob_k.A[ibin][ksum]+=count; 
		if(ibin>=0 && ibin <nbin) prob_k[ibin]+=(count); 

		//ibin=(alph+180.0)/da;
		ibin=(alph+270.0)/da;
		if(ibin>=0 && ibin <nbin) prob_a[ibin]+=(count);


		//k_mean+=(xi*amp);
		//k_sig+=(xi*xi*amp);
		k_mean+=(abs(xi1)*amp);
		k_sig+=(xi1*xi1*amp);
		a_mean+=(alph*amp);
		a_sig+=(alph*alph*amp);
	}
	}

	k_mean/=asum;
	k_sig/=asum;
	a_mean/=asum;
	a_sig/=asum;
	k_sig=sqrt((k_sig-k_mean*k_mean));
	a_sig=sqrt((a_sig-a_mean*a_mean));
	/*
	asum=1.0;
	for(i=0;i<nbin;i++){
		prob_k[i]/=asum;
		prob_a[i]/=asum;
	}
	*/
};

