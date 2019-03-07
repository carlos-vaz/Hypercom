#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define ANNOUNCER_PROC 0

int myrank, rank_2d, mycoord[2], np, dims_procs[2], num_points, dims_pts[2];
int proc_pts[2], proc_size;


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
	 * File metadata in the form: X dimensions, Y dimensions (CARTESIAN FORM)
	 */
	dims_procs[0] = atoi(argv[3]);
	dims_procs[1] = atoi(argv[2]);
	if(myrank == ANNOUNCER_PROC) printf("Will solve %s with %d (x) by %d (y) processes\n", argv[1], dims_procs[0], dims_procs[1]);
	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_read_all(file, &dims_pts, 2, MPI_INT, MPI_STATUS_IGNORE);
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
	MPI_File_set_view(file, 2*sizeof(int)+(mycoord[0]*proc_pts[0]+mycoord[1]*proc_size*dims_procs[0])*sizeof(double), \
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
	double *T = (double*)calloc(proc_size, sizeof(double));
	if(mycoord[1]==0) {
		printf("Process (%d, %d): I have a southern boundary\n", mycoord[0], mycoord[1]);
		// Do nothing (xe^y == 0)
	}
/*	if(mycoord[1]==dims_procs[1]) {
		printf("Process (%d, %d): I have a northern boundary\n", mycoord[0], mycoord[1]);
		for(int x=0; x<dims_pts[0]; x++)
			T[(proc_pts[1]-1)*proc_pts[0]+x] = v[(proc_pts[1]-1)*proc_pts[0]+x]; // copy from 'v'
	}
*/	if(mycoord[0]==0) {
		printf("Process (%d, %d): I have an eastern boundary\n", mycoord[0], mycoord[1]);
		for(int y=0; y<dims_pts[1]; y++)
			T[y*proc_pts[0]] = v[y*proc_pts[0]]; // copy from 'v'
	}
/*	if(mycoord[0]==dims_procs[0]) {
		printf("Process (%d, %d): I have a western boundary\n", mycoord[0], mycoord[1]);
		for(int y=0; y<dims_pts[1]; y++)
			T[(y+1)*proc_pts[0]-1] = v[(y+1)*proc_pts[0]-1];  // copy from 'v'
	}
*/
	MPI_Finalize();
	return 0;
}
