#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int myrank, rank_2d, coord_2d[2], np, dims_procs[2], num_points, dims_pts[2];
int proc_pts[2];


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

	// Check for correct number of arguments
	if(argc!=4) {
		if(myrank==0)
			printf("ERROR: arguments to main: mpirun -np [#procs] main [sourcefile] [# procs X] [# procs Y] \n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}


	/*
	 * Read the file metadata and check for do-ablility
	 */
	dims_procs[0] = atoi(argv[2]);
	dims_procs[1] = atoi(argv[3]);
	printf("Will solve %s with %d (y) by %d (x) processes\n", argv[1], dims_procs[0], dims_procs[1]);
	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_read_all(file, &dims_pts, 2, MPI_INT, MPI_STATUS_IGNORE);
	printf("FILE READ... DIMENSIONS: %d BY %d\n", dims_pts[0], dims_pts[1]);
	if(dims_pts[0]%dims_procs[0]!=0 || dims_pts[1]%dims_procs[1]!=0) {
		if(myrank == 0)
			printf("Number of points must be divisible by number of procs in that dimension\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}


	/*
	 * Create a cartesian topology. This lets us identify processes by their 
	 * coordinates instead of their ranks. 
	 */
	MPI_Comm comm2d;
	int periodic[2] = {0,0};	// What is periodic?
	MPI_Cart_create(MPI_COMM_WORLD, 2,dims_procs,periodic,1,&comm2d);
	MPI_Cart_coords(comm2d,myrank, 2,coord_2d);
	MPI_Cart_rank(comm2d,coord_2d,&rank_2d);
	printf("I am %d: (%d,%d); originally %d\n",rank_2d,coord_2d[0],coord_2d[1],myrank);


	/*
	 * Read the file contents into vector
	 */
	MPI_Datatype vector; 
	proc_pts[0] = dims_pts[0]/dims_procs[0];
	proc_pts[1] = dims_pts[1]/dims_procs[1];
	printf("# blocks per process = %d\t# pts per proc (X) = %d\t# pts per proc (Y) = %d\n",proc_pts[0], proc_pts[1], proc_pts[0]);
	MPI_Type_vector(proc_pts[1], proc_pts[0], proc_pts[1]*np, MPI_DOUBLE, &vector);
	MPI_Type_commit(&vector);


	MPI_Finalize();
	return 0;
}
