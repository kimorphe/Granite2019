#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "heap.h"

//-------------------------------------------
void Heap::init(int n){
	nsize=n;
	ndat=0;
	hp=(double *)malloc(sizeof(double)*nsize);
	indx=(int *)malloc(sizeof(int)*nsize);
	for(int i=0;i<nsize;i++){
		hp[i]=0.0;
		indx[i]=0;
	}
};

void Heap::add(double val,int ID){
	ndat++;
	if(ndat>nsize){
		printf(" Error!!: number of data exceeds array size !\n");
		printf(" ndat=%d nsize=%d\n",ndat,nsize);
		printf(" --> abort process.");
		ndat--;
		exit(-1);
	};
	hp[ndat-1]=val;
	indx[ndat-1]=ID;
	double tmp;
	int c=ndat;
	int p=ndat/2;
	int itmp;
	if(ndat>1){
	while(hp[p-1]>hp[c-1]){
		tmp=hp[p-1];
		hp[p-1]=hp[c-1];
		hp[c-1]=tmp;

		itmp=indx[p-1];
		indx[p-1]=indx[c-1];
		indx[c-1]=itmp;
		c=p;
		p=c/2;
		if(p==0) break;
	};
	}
}
double Heap::pop(int *i0){
	double tmp,dat0;
	dat0=hp[0];
	*i0=indx[0];
	hp[0]=hp[ndat-1];
	indx[0]=indx[ndat-1];
	ndat--;

	int c,c1,c2,p=1;
	c1=p*2;
	int itmp;
	while(c1<=ndat){
		c2=c1+1;
		c=c1;
		if(c2<=ndat){
			if(hp[c2-1]<hp[c1-1]) c=c2;
		};
		if(hp[p-1] < hp[c-1]) break;
		tmp=hp[p-1];
		hp[p-1]=hp[c-1];
		hp[c-1]=tmp;

		itmp=indx[p-1];
		indx[p-1]=indx[c-1];
		indx[c-1]=itmp;
		p=c;
		c1=p*2;
	};
	return(dat0);
};

void Heap::sort(){
	int i,p,c,c1,c2;
	int n=ndat;
	double tmp;
	int itmp;

	for(i=0;i<ndat;i++){
		tmp=hp[0];
		hp[0]=hp[n-1];
		hp[n-1]=tmp;

		itmp=indx[0];
		indx[0]=indx[n-1];
		indx[n-1]=itmp;
		n--;

		p=1;
		c1=p*2;
		while(c1<=n){
			c2=c1+1;
			c=c1;
			if(c2 <= n){
				if(hp[c2-1]<hp[c1-1]) c=c2;
			};
			if(hp[p-1]<hp[c-1]) break;
			tmp=hp[p-1];
			hp[p-1]=hp[c-1];
			hp[c-1]=tmp;

			itmp=indx[p-1];
			indx[p-1]=indx[c-1];
			indx[c-1]=itmp;
			p=c;
			c1=p*2;
		}
	};		
};
void Heap::print(){
	printf("ndat=%d nsize=%d\n",ndat,nsize);
	for(int i=0;i<ndat;i++) printf("%d %d %lf\n",i,indx[i],hp[i]);
	puts("----------------------------");
};
int main(){
	Heap hp;
	int i,nsize=22,ndat=21;
	hp.init(nsize);
	double val;
	for(i=0;i<ndat;i++){
		val=sin(double(i));
		hp.add(val,2*i+1);
	};
	hp.print();

	int i0;
	double dat0;
	puts("pop data one-by-one:");
	for(i=0;i<ndat;i++){
		dat0=hp.pop(&i0);
	       	printf("%lf %d\n",dat0,i0);
	}
	printf("ndat=%d\n",hp.ndat);
	puts("-------------------");

	for(i=0;i<ndat;i++){
		val=sin(double(i));
		hp.add(val,2*i+1);
	};
	hp.print();
	hp.sort();
	hp.print();

	return(0);
};
