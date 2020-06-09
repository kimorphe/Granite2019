#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Heap{
	public:
		//double *dat;	// data 
		double *hp;	// heap 
		int *indx;	// index
		int nsize;	// array size
		int ndat;	// number of data
		void init(int n); 
		void print();
		//void heap();
		void sort();
		double pop(int *i0);
		void  add2heap(double val);
		void add(double val, int ID);
	private:

};
