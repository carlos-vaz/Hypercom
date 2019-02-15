#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX_PROCS 1000

int myrank, np, per_proc, rank_of_children[2] = {-1,-1};

/* Compute single function value. 
 * To swap function, modify this only. 
 */
double func_val(double x) {
	return 4/(1+pow(x,2));
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
 * RETURN VALUE: the number of processes I pro-
 * pagated to (0 if leaf process). 
 */
int propagate(double *arr, int count, double *myshare) {
	// Reached leaf-node process
	if(count == per_proc) {
		memcpy(myshare, arr, per_proc*sizeof(double));
		free(arr);
		return 0;
	}

	// Calculate partitions and destinations
	int chunks = count/per_proc;
	int chunks_left = chunks>>1;
	int chunks_right = chunks-chunks_left-1;
	int dest_rank_left  = myrank - ((chunks_left-1)>>1) - 1;
	int dest_rank_right = myrank + (chunks_right>>1)  + 1;
	int children=0;

	// Allocate left & right sub-arrays
	int larr_sz = chunks_left*per_proc, rarr_sz = chunks_right*per_proc;
	double * larr = malloc(larr_sz*sizeof(double));
	double * rarr = malloc(rarr_sz*sizeof(double));

	// Copy from array to sub-arrays
	memcpy(larr, arr, larr_sz*sizeof(double));
	memcpy(rarr, &arr[larr_sz+per_proc], rarr_sz*sizeof(double));
	memcpy(myshare, &arr[larr_sz], per_proc*sizeof(double));

	// Send work to left and right
	if(chunks_left > 0) {
		children++;
		rank_of_children[0] = dest_rank_left;
		MPI_Send(larr, larr_sz, MPI_DOUBLE, dest_rank_left, 0, MPI_COMM_WORLD);
	}
	if(chunks_right > 0) {
		children++;
		rank_of_children[ rank_of_children[0]==-1 ? 0 : 1 ] = dest_rank_right;
		MPI_Send(rarr, rarr_sz, MPI_DOUBLE, dest_rank_right, 0, MPI_COMM_WORLD);
	}
	
	// Free arrays
	free(arr);
	free(larr);
	free(rarr);
	return children;
}

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	struct timeval stop, start;

	// Parse arguments
	if(argc!=4) {
		printf("Abort. main needs 3 arguments: \"mpirun -np [# procs] main [# Points] [X_min] [X_max]\"\n");
		goto finish;
	}
	int Points = strtol(argv[1], NULL, 10);
	double X_min = atof(argv[2]);
	double X_max = atof(argv[3]);

	// Check divisibility of labor
	if(Points % np != 0) {
		printf("Abort. Number of processes (%d) must divide number of points (%d).\n", np, Points);
		goto finish;
	}
	per_proc = Points/np;

	// Allocate arrays for receiving 
	double * arr;
	int arr_sz=0;
	double * myshare = malloc(per_proc*sizeof(double));
	int children;

	// Remember who sent you data. They will be expecting your computation result
	int rank_of_parent = -1;

	if(myrank == np>>1) {
		// Middle process generates function data, then propagates down tree
		// Middle process starts wall clock timer
		gettimeofday(&start, NULL);
		arr = malloc(Points*sizeof(double));
		func_gen(arr, X_min, X_max, Points);
		children = propagate(arr, Points, myshare);
	} else {
		// Block until you receive a message, then receive and propagate down tree
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &arr_sz);
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &arr_sz);
		arr = malloc(arr_sz*sizeof(double));
		MPI_Recv(arr, arr_sz, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		rank_of_parent = status.MPI_SOURCE;
		children = propagate(arr, arr_sz, myshare);
	}

	// Perform the integration
	double sum=0;
	double Delta = (X_max-X_min)/Points;
	for(int i=0; i<per_proc; i++)
		sum += myshare[i]*Delta;

	// If I am NOT leaf, first receive from my children
	double child_sum=0;
	for(int i=0; i<children; i++) {
		MPI_Recv(&child_sum, 1, MPI_DOUBLE, rank_of_children[i], 0, MPI_COMM_WORLD, &status);
		sum += child_sum;
	}
	
	// Send my sum to my parent (if I am not root)
	if(rank_of_parent != -1)
		MPI_Send(&sum, 1, MPI_DOUBLE, rank_of_parent, 0, MPI_COMM_WORLD);

	printf("PROC %d REPORTS SUM = %lf", myrank, sum);
	if(rank_of_parent ==-1) {
		gettimeofday(&stop, NULL);
		printf("\t<-- Result. ELAPSED TIME: %f sec", (double)(stop.tv_usec-start.tv_usec)/1000000 + stop.tv_sec-start.tv_sec);
	}
	printf("\n");

	free(myshare);

finish:	MPI_Finalize();
	return 0;
}
