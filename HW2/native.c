#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX_PROCS 1000

int myrank, np;
long per_proc;

/* Compute single function value. 
 * To swap function, modify this only. 
 */
double func_val(double x) {
	return 4/(1+pow(x,2));
}

/* Store function values in array. 
 * Only one special process will call
 * this function, and pass the chunks
 * to its friends. 
 */
void func_gen(double * arr, double X_min, double X_max, long Points) {
	double x = X_min, Delta = (X_max-X_min)/Points;
	for(int i=0; i<Points; i++, x += Delta)
		arr[i] = func_val(x);
}

int main(long argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	struct timeval stop, start;

	// Parse arguments
	if(myrank==0 && argc!=4) {
		printf("Abort. main needs 3 arguments: \"mpirun -np [# procs] main [# Points] [X_min] [X_max]\"\n");
		MPI_Abort(MPI_COMM_WORLD, 16);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	long Points = strtol(argv[1], NULL, 10);
	double X_min = atof(argv[2]);
	double X_max = atof(argv[3]);

	// Check divisibility of labor
	if(myrank==0 && (Points % np != 0 || Points == np) ) {
		printf("Abort. Number of processes (%d) must divide number of points (%d), and " \
		"cannot equal number of points.\n", np, Points);
		MPI_Abort(MPI_COMM_WORLD, 17);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	per_proc = Points/np;

	// Propagate data down the tree using MPI natives
	double * arr = NULL;
	double * myshare = malloc(per_proc*sizeof(double));
	if(myrank==0) {
		arr = malloc(Points*sizeof(double));
		func_gen(arr, X_min, X_max, Points);
		gettimeofday(&start, NULL);	// Start timer AFTER you compute function
	}
	printf("Proc0 has already gen %li points\n", Points);
	MPI_Scatter(arr, (int)per_proc, MPI_DOUBLE, myshare, (int)per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	free(arr);

	// Perform the integration
	double sum=0;
	double Delta = (X_max-X_min)/Points;
	for(int i=0; i<per_proc-1; i++)
		sum += (myshare[i]+myshare[i+1])*Delta/2;
	free(myshare);


	// Collect sums up the tree using MPI natives
	double * sums = NULL;
	if(myrank==0)
		sums = malloc(np*sizeof(double));
	MPI_Gather(&sum, 1, MPI_DOUBLE, sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Process 0: collect sums into final sum
	if(myrank==0)
		for(int i=1; i<np; i++)
			sum += sums[i];
	free(sums);

	// Print results
	if(myrank==0) {
		gettimeofday(&stop, NULL);
		printf("PROC %d REPORTS SUM = %lf", myrank, sum);
		printf("\t<-- Result. ELAPSED TIME: %f sec\n", (double)(stop.tv_usec-start.tv_usec)/1000000 + stop.tv_sec-start.tv_sec);
	}

	MPI_Finalize();
	return 0;
}
