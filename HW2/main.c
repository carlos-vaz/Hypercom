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
	if(count == per_proc)
		return;

	// Calculate partitions and destinations
	int chunks = count/per_proc;
	int chunks_left = chunks>>1;
	int chunks_right = chunks-chunks_left-1;
	int dest_rank_left  = myrank + (chunks_left-1)>>1 + 1;
	int dest_rank_right = myrank + (chunks_right>>1)  + 1;

	// Allocate left & right sub-arrays
	int bytes_larr = chunks_left*per_proc, bytes_rarr = chunks_right*per_proc;
	double * larr = malloc(bytes_larr*sizeof(double));
	double * rarr = malloc(bytes_rarr*sizeof(double));

	// Copy from array to sub-arrays
	memcpy(larr, arr, bytes_rarr);
	memcpy(rarr, arr+bytes_larr+per_proc, bytes_rarr);
	memcpy(myshare, arr+bytes_larr, per_proc);

	// Send work to left and right
	MPI_Send(larr, chunks_left*per_proc, MPI_DOUBLE, dest_rank_left, 0, MPI_COMM_WORLD);
	MPI_Send(rarr, chunks_right*per_proc, MPI_DOUBLE, dest_rank_right, 0, MPI_COMM_WORLD);
	
	// Free big array
	free(arr);
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// Parse arguments
	if(argc!=3) {
		printf("Abort. Need 3 arguments: \"mpirun -np [# procs] main [# Points] [X_min] [X_max]\"\n");
		goto finish;
	}
	int Points = *argv[1];
	double X_min = *argv[2];
	double X_max = *argv[3];
	double Delta = (X_max-X_min)/Points;

	// Check divisibility of labor
	if(Points % np != 0) {
		printf("Abort. Number of processes (%d) must divide number of points (%d).\n", np, Points);
		printf("argv[1] = %d\n", *argv[1]);
		printf("argv[1] = %d, argv[2] = %d, argv[3] = %d\n", *argv[1]), *argv[2], *argv[3]);
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
		MPI_Recv(arr, arr_sz, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		propagate(arr, arr_sz, myshare);
	}

	// This process' share of computation
	double sum;
	for(int i=0; i<per_proc; i++)
		sum += myshare[i];
	printf("Process %d got sum = %lf.\n", myrank, sum);
	
	free(myshare);


finish:	MPI_Finalize();
	return 0;
}
