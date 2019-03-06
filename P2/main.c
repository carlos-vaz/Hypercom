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

	
	MPI_Comm comm2d;
	int ndim = 2;
	const int periodic[2] = {0,0};
	//periodic[0] = periodic[1] = 0;				// what is periodic?
	const int dimensions[2] = {5,2};
	//dimensions[0] = 10; dimensions[1] = 10;			// 10 x 10?
	int coord_2d[2];
	int rank_2d;
	MPI_Cart_create(MPI_COMM_WORLD,ndim,dimensions,periodic,1,&comm2d);
	MPI_Cart_coords(comm2d,myrank,ndim,coord_2d);
	MPI_Cart_rank(comm2d,coord_2d,&rank_2d);
	printf("I am %d: (%d,%d); originally %d\n",rank_2d,coord_2d[0],coord_2d[1],myrank);


	return 0;
}
