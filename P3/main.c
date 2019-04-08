#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define XRANGE 1
#define YRANGE 1

long grid_size;
int Px, Py, Tx, Ty;

extern 
void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);


void print_grid(double * G, long size, int sizex) {
	for(int i=0; i<size; i++) {
		printf("%lf, ", G[i]);
		if(i%sizex==sizex-1)
			printf("...\n");
	}
}

int main(int argc, char** argv) {
	if(argc!=5) {
		printf("Usage:\n\t./main [# points X] [# points Y] [# threads X] [# threads Y]\n");
		exit(-1);
	}
	Px = atoi(argv[1]);
	Py = atoi(argv[2]);
	Tx = atoi(argv[3]);
	Ty = atoi(argv[4]);
	grid_size = (long)Px*(long)Py; 
	if((Px%Tx!=0) || (Py%Ty!=0)) {
		printf("Threads in X dimension must divide Points in X dimension, and likewise for Y dimension\n");
		exit(-2);
	}
	printf("Running with parameters:\n\tPoints X = %d, Threads X = %d\n\tPoints Y = %d, Threads Y = %d\n", Px, Tx, Py, Ty);

	// Make Grids
	double * S = (double * )malloc(grid_size*sizeof(double));
	double * T = (double * )calloc(grid_size, sizeof(double));
	for(long i=0; i<grid_size; i++) {
		double val = (i%Px)*((double)XRANGE/Px)*powf(2.718281828, (i/Px)*((double)YRANGE/Py));
		S[i] = val;
		if(i/Px==0 || i/Px==Py-1 || i%Px==0 || i%Px==Px-1)
			T[i] = val;
	}
	//print_grid(T, grid_size, Px);

	print_grid(T, grid_size, Px);

	omp_set_num_threads(Tx*Ty);
	int id, i;
	#pragma omp parallel default(none) shared(S,T,grid_size,Tx,Ty,Px,Py) private(id,i)
	{
		id = omp_get_thread_num();
		printf("Hello from thread %d\n", id);

		#pragma omp for
			for(i=0; i<grid_size; i++) {
				printf("\tSubroutine from thread %d (i=%d, Px = %d)\n", id, i, Px);
				if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1) {
					T[i] = (-1*S[i]*pow(((double)(XRANGE/Px)*(double)(YRANGE/Py)),2)+(T[i-1]+T[i+1])*pow((double)(YRANGE/Py),2)+ \
					(T[i-Px]+T[i+Px])*pow((double)(XRANGE/Px),2))/(2*pow((double)(XRANGE/Px),2)+2*pow((double)(YRANGE/Py),2));

				}
			}
	}

	double Xmin = 0;
	double Ymin = 0;
	double Xmax = XRANGE;	
	double Ymax = YRANGE;
	//VTK_out(Px, Py, &Xmin, &Xmax, &Ymin, &Ymax, T, 0);
	print_grid(T, grid_size, Px);




	free(S);
	free(T);
	return 0;
}
