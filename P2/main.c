#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int myrank, np, x_proc, y_proc, num_points, dimensions_pts[2];

/*void read_file_into_vector(MPI_File *file, void *buf, int count, MPI_Dataype vector, int int_or_double) {
	MPI_File_read_all(*file, buf, count, MPI_
}*/

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(argc!=4) {
		printf("ERROR: arguments to main: mpirun -np [#procs] main [sourcefile] [# procs X] [# procs Y] \n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	x_proc = atoi(argv[2]);
	y_proc = atoi(argv[3]);
	printf("Will solve %s with %d horizontal processes and %d vertical processes\n", argv[1], x_proc, y_proc);

	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_read_all(file, &dimensions_pts, 2, MPI_INT, MPI_STATUS_IGNORE);

	printf("FILE READ... DIMENSIONS: %d BY %d\n", dimensions_pts[0], dimensions_pts[1]);


	MPI_Datatype vector; 
	//MPI_Vector_type();

	MPI_Comm comm2d;
	int ndim = 2;
	int periodic[2] = {0,0};	// What is periodic?
	int dimensions[2];
	dimensions[0] = 1; dimensions[1] = 5;
	int coord_2d[2];
	int rank_2d;
	MPI_Cart_create(MPI_COMM_WORLD,ndim,dimensions,periodic,1,&comm2d);
	MPI_Cart_coords(comm2d,myrank,ndim,coord_2d);
	MPI_Cart_rank(comm2d,coord_2d,&rank_2d);
	printf("I am %d: (%d,%d); originally %d\n",rank_2d,coord_2d[0],coord_2d[1],myrank);


	return 0;
}
