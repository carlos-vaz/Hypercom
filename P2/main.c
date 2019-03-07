#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define ERROR_THRESH  	1e-12	// Max tolerable error (infinity norm)
#define ANNOUNCER_PROC 	0
#define	   index(x,y)	x+y*proc_pts[0]
#define	   up(i)	i+proc_pts[0]
#define	   down(i)	i-proc_pts[0]
#define    right(i)	i+1;
#define    left(i)	i-1;

int myrank, rank_2d, mycoord[2], np, dims_procs[2], num_points, dims_pts[2], proc_pts[2], proc_size \
	rank_right, rank_left, rank_up, rank_down;
double deltas[2];


/*
 * Each process is assigned a grid section via MPI_Cart_create, and the corresponding 
 * region of the data file is read into each process using MPI_Type_vector and 
 * MPI_File_read_all. 
 *  
 */
int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	/*
	 * Check for correct number of arguments
	 * Arguments in the form: Y procs, X procs (MATRIX FORM)
	 */
	if(argc!=4) {
		if(myrank == ANNOUNCER_PROC)
			printf("ERROR: arguments to main: mpirun -np [#procs] main [sourcefile] [# procs Y] [# procs X] \n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}


	/*
	 * Read the file metadata and check for do-ablility
	 * File metadata in the form: X dimensions, Y dimensions, X range, Y range (CARTESIAN FORM)
	 */
	dims_procs[0] = atoi(argv[3]);
	dims_procs[1] = atoi(argv[2]);
	if(myrank == ANNOUNCER_PROC) printf("Will solve %s with %d (x) by %d (y) processes\n", argv[1], dims_procs[0], dims_procs[1]);
	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_read_all(file, &dims_pts, 2, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read_all_at(file, 2*sizeof(int), &deltas, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
	if(myrank == ANNOUNCER_PROC) printf("FILE READ... DIMENSIONS (MATRIX FORM): %d BY %d\n", dims_pts[1], dims_pts[0]);
	if( (dims_pts[0]%dims_procs[0]!=0 || dims_pts[1]%dims_procs[1]!=0) && myrank==ANNOUNCER_PROC) {
		printf("Number of points must be divisible by number of procs in that dimension.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	if(np!=dims_procs[0]*dims_procs[1] && myrank == ANNOUNCER_PROC) {
		printf("Allocated %d processes, but user specified %d*%d = %d processes.\n", \
				np, dims_procs[1], dims_procs[0], dims_procs[1]*dims_procs[0]);
		MPI_Abort(MPI_COMM_WORLD, -1);	
	}
	MPI_Barrier(MPI_COMM_WORLD); // wait in case proc 0 aborts


	/*
	 * Create a cartesian topology. This lets us identify processes by their 
	 * coordinates instead of their ranks. 
	 */
	MPI_Comm comm2d;
	int periodic[2] = {0,0};	// What is periodic?
	MPI_Cart_create(MPI_COMM_WORLD, 2,dims_procs,periodic,1,&comm2d);
	MPI_Cart_coords(comm2d,myrank, 2,mycoord);
	MPI_Cart_rank(comm2d,mycoord,&rank_2d);
	printf("I am %d: (%d,%d); originally %d\n",rank_2d,mycoord[0],mycoord[1],myrank);
	


	/*
	 * Read the file contents into vector
	 */
	MPI_Datatype vector; 
	proc_pts[0] = dims_pts[0]/dims_procs[0]; // # pts in Y dim of process partition
	proc_pts[1] = dims_pts[1]/dims_procs[1]; // # pts in X dim of process partition
	proc_size = proc_pts[0]*proc_pts[1];
	double *v = (double*)malloc(proc_size*sizeof(double));
	if( myrank == ANNOUNCER_PROC) printf("# blocks per process = %d\t# pts per proc (X) = %d\t# pts per proc (Y) = %d\n", \
							proc_pts[1], proc_pts[0], proc_pts[1]);
	MPI_Type_vector(proc_pts[1], proc_pts[0], dims_pts[0], MPI_DOUBLE, &vector);
	MPI_Type_commit(&vector);	
	printf("MPI_File_set_view(file, 2*sizeof(int)+%d*%d+%d*%d, MPI_DOUBLE, vector, NULL, MPI_INFO_NULL)\n", \
								mycoord[0], proc_pts[0] ,mycoord[1], proc_size);
	MPI_File_set_view(file, 2*sizeof(int)+2*sizeof(double)+(mycoord[0]*proc_pts[0]+mycoord[1]*proc_size*dims_procs[0])*sizeof(double), \
								MPI_DOUBLE, vector, "native", MPI_INFO_NULL);	
	MPI_File_read_all(file, v, proc_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0; i<np; i++) {
		if(myrank==i) {
			printf("\nRANK %d (%d, %d)\n", myrank, mycoord[0], mycoord[1]);
			for(int j=0; j<proc_size; j++) {
				printf("%d, ", (int)v[j]);
			}
			printf("\n\n");
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}


	/*
	 * Each process creates a temperature vector T and fills it with 
	 * 0, except for at the boundaries, which take the initial value 
	 * of xe^y (so we copy from vector 'v' at boundary points)
	 */
	int bound_south = 0;
	int bound_north = 0;
	int bound_east = 0;
	int bound_west = 0;
	double *T = (double*)calloc(proc_size, sizeof(double));
	if(mycoord[1]==0) {
		printf("Process (%d, %d): I have a southern boundary\n", mycoord[0], mycoord[1]);
		bound_south = 1;
		// Do nothing (xe^y == 0)
	}
	if(mycoord[1]==dims_procs[1]-1) {
		printf("Process (%d, %d): I have a northern boundary\n", mycoord[0], mycoord[1]);
		bound_north = 1;
		for(int x=0; x<proc_pts[0]; x++)
			T[(proc_pts[1]-1)*proc_pts[0]+x] = v[(proc_pts[1]-1)*proc_pts[0]+x]; // copy from 'v'
	}
	if(mycoord[0]==0) {
		printf("Process (%d, %d): I have an eastern boundary\n", mycoord[0], mycoord[1]);
		bound_east = 1;
		for(int y=0; y<proc_pts[1]; y++)
			T[y*proc_pts[0]] = v[y*proc_pts[0]]; // copy from 'v'
	}
	if(mycoord[0]==dims_procs[0]-1) {
		printf("Process (%d, %d): I have a western boundary\n", mycoord[0], mycoord[1]);
		bound_west = 1;
		for(int y=0; y<proc_pts[1]; y++)
			T[(y+1)*proc_pts[0]-1] = v[(y+1)*proc_pts[0]-1];  // copy from 'v'
	}


	// Preparation complete... Start the timer

	double *send_south = (double*)malloc(proc_pts[0]*sizeof(double));
	double *send_north = (double*)malloc(proc_pts[0]*sizeof(double));
	double *send_east  = (double*)malloc(proc_pts[1]*sizeof(double));
	double *send_west  = (double*)malloc(proc_pts[1]*sizeof(double));
	double *recv_south = (double*)malloc(proc_pts[0]*sizeof(double));
	double *recv_north = (double*)malloc(proc_pts[0]*sizeof(double));
	double *recv_east  = (double*)malloc(proc_pts[1]*sizeof(double));
	double *recv_west  = (double*)malloc(proc_pts[1]*sizeof(double));
	int got_south = 0;
	int got_north = 0;
	int got_east  = 0;
	int got_west  = 0;
	double myError = ERROR_THRESH + 1;

	MPI_Cart_shift(comm2d, 0, +1, &rank_2d, &rank_right);
	MPI_Cart_shift(comm2d, 0, -1, &rank_2d, &rank_left);
	MPI_Cart_shift(comm2d, 1, +1, &rank_2d, &rank_up);
	MPI_Cart_shift(comm2d, 1, -1, &rank_2d, &rank_down);
	printf("RANK %d (%d,%d): rank_right=%d, rank_up=%d\n", myrank, mycoords[0], mycoords[1], rank_right, rank_left);

//	while(myError > ERROR_THRESH) {
		/*
		 * Post a non-blocking send and a non-blocking receive to all neighbors.
		 * While you update your internal temperatures, hopefully the requests
		 * will go through. 
		 */
/*		if(!border_south) {
			MPI_Irecv(recv_south, proc_pts[0], MPI_DOUBLE, rank_down, );
			MPI_Isend();
		}
*/
		/*
		 * Compute & update temperatures on block interiors (Gauss-Seidel) for
		 * 1 iteration. 
		 */
/*		int i=0;
		for(int y=1; y<dims_procs[1]-1; y++) {
			for(int x=1; x<dims_procs[0]-1; x++) {
				i = index(x,y);
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+T[right(i)])*pow(deltas[1],2)+ \
					(T[down(i)]+T[up(i)]*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
*/		}

		/*
		 * Check the status of your send and receive requests. 
		 */

//	}


	free(v);
	free(T);
	free(send_south);
	free(send_north);
	free(send_east);
	free(send_west);
	free(recv_south);
	free(recv_north);
	free(recv_east);
	free(recv_west);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
