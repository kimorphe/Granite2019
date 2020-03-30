#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#ifndef __WVS__
	#define __WVS__
	#include "waves.h"
#endif

using namespace std;

//------------------------------------------------------------
Wv1D::Wv1D(){
	//amp=0;
	//time=0;
	mllc=false;
	fft_stat=0;
};
Wv1D::Wv1D(int ndat){
	Nt=ndat;
	amp=(double *)malloc(sizeof(double)*Nt);
	time=(double *)malloc(sizeof(double)*Nt);
	mllc=true;
	fft_stat=0;
	strcpy(data_file,"");
	t1=0.0;
	t2=0.0;
	dt=0.0;
};
Wv1D::Wv1D(char *fname){
	amp=0;
	time=0;
	mllc=false;
	fft_stat=0;
	Wv1D::load(fname);
};
int Wv1D::count_lines(char *fname){
	FILE *fp=fopen(fname,"r");
	char cbff[128];
	Nt=0;
	while(fgets(cbff,128,fp)!=NULL){
		Nt++;
	}
	fclose(fp);
	return(Nt);
};
int Wv1D::load2(char *fname,bool count){

	int nn;
	if(count) nn=Wv1D::count_lines(fname);

	FILE *fp=fopen(fname,"r");
	strcpy(data_file,fname);
	if(!mllc){
		amp=(double *)malloc(sizeof(double)*Nt);
		mllc=true;
	}
	time=(double *)malloc(sizeof(double)*Nt);
	double sum=0.0;
	for(int i=0;i<Nt;i++){
		fscanf(fp,"%lf, %lf\n",time+i,amp+i);
		time[i]*=1.e06;
		sum+=amp[i];
	};
	t1=time[0];
	t2=time[Nt-1];
	dt=time[1]-time[0];
	sum/=Nt;
	for(int i=0;i<Nt;i++) amp[i]-=sum;
	//printf("Nt=%d\n",nn);
	//printf("sum=%lf\n",sum);
	fclose(fp);
};
void Wv1D::print_info(){
	printf("File: %s\n",data_file);
	printf("(t1,t2)=%lf,%lf\n",t1,t2);
	printf("dt=%lf\n",dt);
	printf("Nt=%d\n",Nt);
}
int Wv1D::load(char *fname){
	FILE *fp=fopen(fname,"r");

	strcpy(data_file,fname);
	fscanf(fp,"%lf, %lf\n",&t1,&dt);
	t1*=1.e06;
	dt*=1.e06;
	fscanf(fp,"%d\n",&Nt);
	if(!mllc){
		amp=(double *)malloc(sizeof(double)*Nt);
		time=(double *)malloc(sizeof(double)*Nt);
		mllc=true;
	}
	double sum=0.0;
	for(int i=0;i<Nt;i++){
		fscanf(fp,"%lf\n",amp+i);
		time[i]=t1+dt*i;
		sum+=amp[i];
	};
	t2=time[Nt-1];
	sum/=Nt;
	for(int i=0;i<Nt;i++) amp[i]-=sum;
	fclose(fp);
}
int Wv1D::FFT(int isgn){
	int p=ceil(log2(Nt));
	Np=pow(2,p);

	if(fft_stat==0) Amp=(complex<double> *)malloc(sizeof(complex<double>)*Np);
	if(isgn==1){
		for(int i=0;i<Nt;i++) Amp[i]=complex<double>(amp[i],0.0);
		for(int i=Nt;i<Np;i++) Amp[i]=complex<double>(0.0,0.0);
	}

	fft(Amp,Np,isgn);
	fft_stat=isgn;
	df=1./dt/Np;

	return(Np);
};
double Wv1D::L2(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double sum=0.0;
	for(i=i1;i<=i2;i++) sum+=(amp[i]*amp[i]);
	return(sqrt(sum)/(i2-i1+1));
};
double Wv1D::max(double t1, double t2){
	int i1,i2;
	i1=int(t1/dt);
	i2=int(t2/dt);
	if(i1<0) i1=0;
	if(i2<0) i2=0;
	if(i1>=Nt-1) i1=Nt-1;
	if(i2>=Nt-1) i2=Nt-1;

	int i=0;
	double amax=fabs(amp[i2]);
	for(i=i1;i<i2;i++){
	       	if(amax < fabs(amp[i])) amax=fabs(amp[i]);
	}
	return(amax);
};
void Wv1D::out_Amp(char *fn,int ofst){
	FILE *fp=fopen(fn,"w");
	int j;
	double xx,x1,dx;
	if(fft_stat==1){
		dx=df;
		x1=0.0;
	}else{
		x1=t1;
		dx=dt;
	}
	for(int i=0;i<Np;i++){
		j=(ofst+i)%Np;
		xx=(i-ofst)*dx;
		fprintf(fp,"%le %le %le %le\n",xx,Amp[j].real(),Amp[j].imag(),abs(Amp[j]));
	}
	fclose(fp);
};
void Wv1D::out_amp(char *fn){
	FILE *fp=fopen(fn,"w");
	for(int i=0;i<Nt;i++) fprintf(fp,"%le %le\n",i*dt+t1,amp[i]);

	fclose(fp);
};
void Wv1D::Butterworth(double tb, double Tw_6dB){

	int p=4;
	double tt;
	double t0=Tw_6dB*0.5;
	double arg;
	for(int i=0; i<Nt;i++){
		tt=t1+dt*i;
		arg=(tt-tb)/t0;
		arg=pow(arg,p);
		amp[i]/=(1.+arg);
	};
};
double Wv1D::gdelay(){
	if(fft_stat==0) FFT(1);
	tg=(double *)malloc(sizeof(double)*(Nt-1));
	double *phi=(double *)malloc(sizeof(double)*Nt);
	double *tmp=(double *)malloc(sizeof(double)*Nt);

	int i;
	double phi0=arg(Amp[0]),phi1;
	double PI2=8.0*atan(1.0);
	double PI=4.0*atan(1.0);

	// phase spectrum
	double cost,sint,zlen;
	double Amax=0.0;
	for(i=0;i<Nt;i++){
		zlen=abs(Amp[i]);
		cost=real(Amp[i])/zlen;	
		sint=imag(Amp[i])/zlen;	
		tmp[i]=acos(cost);
		if( sint < 0.0) tmp[i]=PI2-tmp[i];
		if(Amax<abs(Amp[i])) Amax=abs(Amp[i]);
	};

	// unwrapping
	double dphi;
	phi[0]=tmp[0];
	for(i=1;i<Nt;i++){
		dphi=tmp[i]-tmp[i-1];
		if(dphi>PI) dphi=-(PI2-dphi);
		if(dphi<-PI) dphi=(PI2+dphi);
		phi[i]=phi[i-1]+dphi;
	};

	double eps=1.e-02;
	double Wf;	//Weiner Filter	
	double Ps,Pn;
	Pn=Amax*eps;
	Pn*=Pn;
	for(i=1;i<Nt-1;i++){
		Ps=abs(Amp[i]*conj(Amp[i]));
		Wf=Ps/(Ps+Pn);
		Amp[i]*=Wf;
	}

	// group delay 
	double df=1./dt/Np;
	double Cf=PI2*df;
	double f1=0.05, f2=1.0;
	FILE *fp=fopen("tg.out","w");
	int nf1=int(f1/df), nf2=int(f2/df);
	for(i=0;i<Nt-1;i++){
		tg[i]=(phi[i+1]-phi[i])/Cf;
		fprintf(fp,"%lf %lf %lf %lf\n",df*i,tg[i],1./tg[i],phi[i]);
	};
	fclose(fp);
	double tgb=0.0;
	for(i=nf1;i<=nf2;i++) tgb+=tg[i];
	return(tgb/(nf2-nf1+1));
};
//------------------------------------------------------------
Wv1D corr(Wv1D wv1, Wv1D wv2, double *tmax, double *Amax){
	wv1.FFT(1);
	wv2.FFT(1);

	Wv1D wv3(wv1.Np);
	wv3.FFT(1);
	double A1=0.0,A2=0.0;
	complex <double > Z1,Z2;
	Z1=complex<double>(0.0,0.0);
	Z2=complex<double>(0.0,0.0);
	wv3=wv1;
	for(int i=0;i<wv1.Np;i++){
		Z1=wv1.Amp[i];
		Z2=conj(wv2.Amp[i]);
		wv3.Amp[i]=Z1*Z2;
		A1+=abs(Z1*Z1);
		A2+=abs(Z2*Z2);
	};
	wv3.FFT(-1);
	A1=sqrt(A1*A2);
	int imax=0;
	double A0=abs(wv3.Amp[0].real()/A1);
	for(int i=0;i<wv1.Np;i++){
		wv3.Amp[i]/=A1;
		if(A0 < abs(wv3.Amp[i].real())){
			A0=abs(wv3.Amp[i].real());
			imax=i;
		}
	}
	double t0=imax*wv1.dt;
	if(imax>wv1.Np/2) t0=-(wv1.Np-imax)*wv1.dt;
	(*tmax)=t0;
	(*Amax)=A0;
	//printf("%d %le %le\n",imax,t0,A0);
	return(wv3);
};
//---------------------------------------------------------
#if DB == 11 
int main(){
	char fname[128];
	sprintf(fname,"../W20H30_fine/scope_%d.csv",80);
	Wv1D awv(fname);
	//awv.Butterworth(16.0,10.0);
	awv.FFT(1);

	// Output
	sprintf(fname,"awvt.out");
	awv.out_amp(fname);
	sprintf(fname,"awvw.out");
	awv.out_Amp(fname,0);
	printf("tg=%lf\n",awv.gdelay());
	return(0);
};
#endif
