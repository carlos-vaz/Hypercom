#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX_PROCS	1000

int myrank, np, per_proc;

/* Compute single function value. 
 * To swap function, modify this only. 
 */
double func_val(double x) {
	return x;
}

/* Store function values in array. 
 * Only one special process will call
 * this function, and pass the chunks
 * to its friends
 */
void func_gen(double * arr, double X_min, double X_max, int Points) {
	double x = X_min;
	double Delta = (X_max-X_min)/Points;
	for(int i=0; i<Points; i++) {
		arr[i] = func_val(x);
		x += Delta;
	}
}

/* Takes as input an array of size count and 
 * an empty array with size per_proc. Will copy
 * its share into myshare, and free memory of
 * arr. It will also propagate the remaining 
 * data to the left & right of its share to 
 * the appropriate process down the binary tree. 
 */
void propagate(double *arr, int count, double *myshare) {
	// Reached leaf-node process
	if(count == per_proc) {
		printf("leaf node (proc %d)\n", myrank);
		memcpy(myshare, arr, per_proc*sizeof(double));
		printf("myshare[%d] = %lfn\n", per_proc-10, myshare[per_proc-10]);
		free(arr);
		return;
	}

	// Calculate partitions and destinations
	int chunks = count/per_proc;
	int chunks_left = chunks>>1;
	int chunks_right = chunks-chunks_left-1;
	int dest_rank_left  = myrank - ((chunks_left-1)>>1) - 1;
	printf("myrank=%d, chunks_left=%d, (chunks_left-1)>>1=%d, dest_left=%d\n", myrank, chunks_left, (chunks_left-1)>>1, dest_rank_left);
	int dest_rank_right = myrank + (chunks_right>>1)  + 1;

	// Allocate left & right sub-arrays
	int larr_sz = chunks_left*per_proc, rarr_sz = chunks_right*per_proc;
	double * larr = malloc(larr_sz*sizeof(double));
	double * rarr = malloc(rarr_sz*sizeof(double));

	// Copy from array to sub-arrays
	printf("arr[490]=%lf\n", arr[490]);
	memcpy(larr, arr, larr_sz*sizeof(double));
	memcpy(rarr, &arr[larr_sz+per_proc], rarr_sz*sizeof(double));
	memcpy(myshare, &arr[larr_sz], per_proc*sizeof(double));

	// Send work to left and right
	printf("(%d) sending to left (%d, %d chunks) and right (%d, %d chunks)\n", myrank, dest_rank_left, chunks_left, dest_rank_right, chunks_right);
	if(chunks_left > 0) {
		printf("(%d) sending larr. larr_sz = %d, larr[%d] = %lf", myrank, larr_sz, per_proc-10, larr[per_proc-10]);
		MPI_Send(larr, larr_sz, MPI_DOUBLE, dest_rank_left, 0, MPI_COMM_WORLD);
	}
	if(chunks_right > 0)
		MPI_Send(rarr, rarr_sz, MPI_DOUBLE, dest_rank_right, 0, MPI_COMM_WORLD);
	
	// Free arrays
	free(arr);
	free(larr);
	free(rarr);
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Parse arguments
	if(argc!=4) {
		printf("Abort. main needs 3 arguments: \"mpirun -np [# procs] main [# Points] [X_min] [X_max]\"\n");
		//printf("argv[1] = %d\n", *argv[1]);
		printf("argv[1] = %d, argv[2] = %lf, argv[3] = %lf\n", (int)strtol(argv[1], NULL, 10), atof(argv[2]), atof(argv[3]));
		goto finish;
	}
	int Points = strtol(argv[1], NULL, 10);
	double X_min = atof(argv[2]);
	double X_max = atof(argv[3]);
	double Delta = (X_max-X_min)/Points;

	// Check divisibility of labor
	if(Points % np != 0) {
		printf("Abort. Number of processes (%d) must divide number of points (%d).\n", np, Points);
		goto finish;
	}
	per_proc = Points/np;

	printf("Process %d reporting for duty.\n", myrank);

	double * arr;
	int arr_sz=0;
	double * myshare = malloc(per_proc*sizeof(double));

	if(myrank == np>>1) {
		// Middle process generates function data
		arr = malloc(Points*sizeof(double));
		func_gen(arr, X_min, X_max, Points);
		propagate(arr, Points, myshare);
	} else {
		// Block until you receive a message
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &arr_sz);
		arr = malloc(arr_sz*sizeof(double));
		printf("(%d) receiving %d chunks\n", myrank, arr_sz/per_proc);
		MPI_Recv(arr, arr_sz, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		propagate(arr, arr_sz, myshare);
	}

	// This process' share of computation
	double sum=0;
	for(int i=0; i<per_proc; i++)
		sum += myshare[i];
	printf("Process %d got sum = %lf.\n", myrank, sum);
	
	free(myshare);


finish:	MPI_Finalize();
	return 0;
}
