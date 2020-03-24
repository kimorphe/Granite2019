#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "waves.h"

using namespace std;

int main(int argc, char *argv[]){
	char dir_name[128];
	char fname[128];
	char fnout[128]="Amp.dat";
	int i,num;
	Wv1D awv;
	sprintf(dir_name,"%s","../CoreM_short3/x30y20");

	for(i=0;i<1;i++){
		sprintf(fname,"%s/scope_%d.csv",dir_name,i);
		puts(fname);
		printf("Nt=%d\n",awv.count_lines(fname));
		awv.load2(fname);
		awv.FFT(1);
		awv.out_Amp(fnout,0);
		awv.gdelay();
	}
	
	return(0);
};
