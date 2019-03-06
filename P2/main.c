#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int myrank, np, num_points;

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(argc!=2) {
		printf("ERROR: arguments to main: mpirun -np [#procs] main [sourcefile]\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, "data.txt", MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_FILE_read_all(file, &num_points, 1, MPI_INT, MPI_STATUS_IGNORE);

	printf("FILE CONTAINS %d ELEMENTS\n", num_points);

	MPI_Datatype vector; 
	MPI_Vector_type();
	
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
