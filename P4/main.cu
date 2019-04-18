#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int xdim, ydim;

int main(int argc, char **argv) {
	if(argc!=3) {
		printf("Usage:\n\t./main [xdim] [ydim]\n");
		exit(0);
	}
	xdim = atoi(argv[1]);
	ydim = atoi(argv[2]);

	// Make Grids
	double * S = (double * )malloc(grid_size*sizeof(double));
	double * T = (double * )calloc(grid_size, sizeof(double));
	double * T_tmp = (double * )calloc(grid_size, sizeof(double));
	for(long i=0; i<grid_size; i++) {
		double val = (i%Px)*((double)XRANGE/Px)*powf(2.718281828, (i/Px)*((double)YRANGE/Py));
		//printf("TEST: %ld/%d = %lf\n", i, Px, (double)(i/Px));
		S[i] = val;
		if(i/Px==0 || i/Px==Py-1 || i%Px==0 || i%Px==Px-1) {
			T[i] = val;
			T_tmp[i] = val;
		}
	}
	
}
