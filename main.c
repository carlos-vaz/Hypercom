#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define POINTS 	1000
#define X_MIN	0.0
#define X_MAX	1.0
#define DELTA	(X_MIN-X_MAX)/POINTS


/* Compute single function value. 
 * To swap function, modify this only. 
 */
double func_val(double x) {
	return pow(x,2);
}

/* Store function values in array. 
 * Only one special process will call
 * this function, and pass the chunks
 * to its friends
 */
void func_gen(double * arr) {
	double x = X_MIN;
	for(int i=0; i<POINTS; i++) {
		arr[i] = func_val(x);
		x += DELTA:
	}
}

int main(int argc, char* argv[]) {
	int myrank, np;
	MPI_Init(&argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	printf("Process %d reporting for duty.\n", myrank);

	double *arr = malloc(POINTS*sizeof(double));
	

	MPI_Finalize();
	return 0;
}
