#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

int nt, np;

double function(double x) {
	return 4/(1+pow(x,2));
}

int main(int argc, char *argv[]) {
	/* 
	 * Parse arguments
	 */
	if(argc!=3) {
		printf("Usage:\n\t./main [# points] [# threads]\n");
		exit(0);
	}
	np = atoi(argv[1]);
	nt = atoi(argv[2]);
	if(np%nt!=0) {
		printf("Error: # threads must divide # points\n");
		exit(1);
	}
	printf("Running %d threads on %d points\n", nt, np);

	/* 
	 * Build function
	 */
	double *buf = malloc(np*sizeof(double));
	double x;
	for(int i=0; i<np; i++) {
		buf[i] = function( (double)i/np );
	}
	
	time_t start = clock();



	return 0;
}
