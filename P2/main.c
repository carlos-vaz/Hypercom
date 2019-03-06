#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int myrank, np;

int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(argc!=3) {
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	
	MPI_Comm comm2d;
	int ndim = 2;
	int periodic[2];
	periodic[0] = 0; periodic[1] = 0;			// what is periodic?
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
