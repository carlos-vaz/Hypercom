#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MAX_PROCS 1000

int myrank, np, per_proc, num_children=0;

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
void func_gen(double * arr, double X_min, double X_max, int Points) {
	double x = X_min, Delta = (X_max-X_min)/Points;
	for(int i=0; i<Points; i++, x += Delta)
		arr[i] = func_val(x);
}

/* Takes as input an array of size count and 
 * an empty array myshare with size per_proc. 
 * Will copy its share into myshare (if this
 * is a leaf-reaching propagate), and free 
 * memory of arr. It will also propagate the
 * remaining data to the right of its share to 
 * the appropriate process down the binary tree.
 * Proactively halves count for next propagation.  
 * 
 * RETURN VALUE: 0 if finished reached leaf node
 * process. 1 if further propagations needed. 
 */
int propagate(double *arr, int *count, double *myshare) {
	// Reached leaf-node process. copy data to myshare
	if(*count == per_proc) {
		memcpy(myshare, arr, per_proc*sizeof(double));
		free(arr);
		return 0;
	}

	// Calculate partitions and destinations
	int chunks = *count/per_proc;
	int send_chunks = chunks>>1;
	int dest = myrank + send_chunks;

	// Allocate send array
	int send_sz = per_proc*send_chunks;
	double * send_arr = malloc(send_sz*sizeof(double));

	// Copy from array to send arrays
	memcpy(send_arr, &arr[*count-send_sz], send_sz*sizeof(double));

	// Send work to next node down the tree
	MPI_Send(send_arr, send_sz, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD);
	
	// Free send_array (arr still needed for next propagate)
	free(send_arr);

	// This proc will expect a computation from each of its children
	num_children++;
	*count >>= 1;
	return 1;
}

int main(int argc, char* argv[]) {
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
	int Points = strtol(argv[1], NULL, 10);
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

	// Allocate arrays for receiving 
	double * arr;
	int arr_sz=0;
	double * myshare = malloc(per_proc*sizeof(double));

	// Remember who sent you data. They will be expecting your computation result
	int rank_of_parent = -1;

	if(myrank == 0) {
		// Middle process generates function data, then propagates down tree
		// Middle process starts wall clock timer
		gettimeofday(&start, NULL);
		arr = malloc(Points*sizeof(double));
		func_gen(arr, X_min, X_max, Points);
		arr_sz = Points;
		while(propagate(arr, &arr_sz, myshare)==1);
	} else {
		// Block until you receive a message, then receive and propagate down tree
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &arr_sz);
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &arr_sz);
		arr = malloc(arr_sz*sizeof(double));
		MPI_Recv(arr, arr_sz, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		rank_of_parent = status.MPI_SOURCE;
		while(propagate(arr, &arr_sz, myshare)==1);
	}

	// Perform the integration
	printf("(%d) INTEGRATING...\n", myrank);
	double sum=0;
	double Delta = (X_max-X_min)/Points;
	for(int i=0; i<per_proc-1; i++)
		sum += (myshare[i]+myshare[i+1])*Delta/2;

	/* //not needed
	double * sums = NULL;
	if(myrank==0)
		sums = malloc(np*sizeof(double));
	MPI_Gather(&sum, 1, MPI_DOUBLE, sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(myrank==0) {
		for(int i=1; i<np; i++)
			sum += sums[i];
		printf("PROC %d: SUM = %lf\n", myrank, sum);
	}
	*/
	MPI_Barrier(MPI_COMM_WORLD);

	// Gather the sums into proc 0
	double recv_sum=0;
	int rank_decay = myrank;
	for(int i, np_grow=2; np_grow<np; np_grow<<=1, rank_decay>>=1, i++) {
		printf("(%d) ITERATION (%d): rank_decay = %d\n", myrank, i, rank_decay);
		if(rank_decay % np_grow == 1) {
			printf("(%d) ITERATION (%d): snd to = %d\n", myrank, i, myrank-(np_grow>>1));
			MPI_Send(&sum, 1, MPI_DOUBLE, myrank-(np_grow>>1), 5, MPI_COMM_WORLD);
			break; // Whoever you sent to now represents your group. you leave. 
		}
		else if(rank_decay % np_grow == 0) {
			printf("(%d) ITERATION (%d): rcv from = %d\n", myrank, i, myrank+(np_grow>>1));
			MPI_Recv(&recv_sum, 1, MPI_DOUBLE, myrank+(np_grow>>1), 5, MPI_COMM_WORLD, &status);
			sum += recv_sum;
		}
	}

	// Display results
	if(myrank==0) {
		gettimeofday(&stop, NULL);
		printf("\nPROC %d REPORTS SUM = %lf\t <-- Result. ELAPSED TIME: %f sec\n", myrank, sum, (double)(stop.tv_usec-start.tv_usec)/1000000 + stop.tv_sec-start.tv_sec);
	}

	free(myshare);

	MPI_Finalize();
	return 0;
}
