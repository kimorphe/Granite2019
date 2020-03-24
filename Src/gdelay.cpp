#define DB 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex>
#include <math.h>
#include "waves.h"

int main(int argc, char *argv[]){
	char dir_name[128];
	char fname[128];
	int i,num;
	Wv1D awv();
	sprintf(dir_name,"%s","../CoreM_short3/x30y20");

	for(i=0;i<2;i++){
		sprintf(fname,"%s/scope_%d.csv",dir_name,i);
		puts(fname);
		awv.load(fname);
	}
	
	return(0);
};
