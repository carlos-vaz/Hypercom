#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define XRANGE 1
#define YRANGE 1

int Px, Py, Tx, Ty, Ptx, Pty, shared_pts;
long grid_size;


extern 
void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);
/*
 * If ran with the index of a shared point (thread-to-thread boundary point) as input, output is 
 * index of the corresponding mutex. 
 */
/*int mutex_map(int i) {
	int map;
	int lX = i%Ptx ? ((i%Px)/Ptx)-1 : (((i+1)%Px)/Ptx)-1;
	int lY = (i/Px)%Pty ? ((i/Px)/Pty)-1 : (((i/Px)+1)/Pty)-1;

	if(i%Ptx==0 || i%Ptx==Ptx-1) {
		// Lies on a Vertical boudnary
		map = lX*(Py-2-lY) + i/Px - lY;
		return map;
	}
	if((i/Px)%Pty==0 || (i/Px)%Pty==Pty-1) {
		// Lies on a Horizontal boudnary
		map = (Tx-1)*(Py-2) lX*(Py-2-lY) + i/Px;
		return map;

	}	
}*/

/*int mutex_map(int i, int id) {
	int per_thread = 2*Ptx+2*Pty;
	int map = per_thread*id;
	int xt = i%Ptx;
	int yt = (i/Px)%Pty;
	if(xt==0 && yt==0) {
		return map + 0;
	}
	if(xt==0 && yt==Pty-1) {
		return map + 0;
	}
	
}*/

int mutex_map(int i, int id) {
	int ret = 1;
	if(i%Ptx==0)
		ret *= (id - 1);
	if(i%Ptx==Ptx-1)
		ret *= (id + 1);
	if((i/Px)%Pty==0)
		ret *= (id - Tx);
	if((i/Px)%Pty==0)
		ret *= (id + Tx);
	return ret;
}

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
	Ptx = Px/Tx;
	Pty = Py/Ty;
	shared_pts = (Tx-1)*(Py-2) + (Ty-1)*(Px-2) - 4*(Tx-1)*(Ty-1);
      //omp_lock_t * mutexes = malloc(shared_pts*sizeof(omp_lock_t)); Lock array of minimum size, but hard to code
	omp_lock_t * mutexes = malloc(Tx*Ty*sizeof(omp_lock_t));
	for(int i=0; i<shared_pts; i++) omp_init_lock(&mutexes[i]);
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
		printf("TEST: %ld/%d = %lf\n", i, Px, (double)(i/Px));
		S[i] = val;
		if(i/Px==0 || i/Px==Py-1 || i%Px==0 || i%Px==Px-1)
			T[i] = val;
	}

	print_grid(T, grid_size, Px);

	omp_set_num_threads(Tx*Ty);
	int id, i, unlock=0;
	#pragma omp parallel default(none) shared(S,T,grid_size,Tx,Ty,Px,Py,Ptx,Pty,mutexes) private(id,i,unlock)
	{
		id = omp_get_thread_num();
		printf("Hello from thread %d\n", id);

		#pragma omp for
			for(i=0; i<grid_size; i++) {
				printf("\tSubroutine from thread %d (i=%d, Px = %d)\n", id, i, Px);
				if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1) {
					if(i%Ptx==Ptx-1 || i%Ptx==0 || (i/Px)%Pty==Pty-1 || (i/Px)%Pty==0) {
						// Updating a shared point. Must lock!
						omp_set_lock(&mutexes[mutex_map(i, id)]);
						unlock = 1;
					}
					T[i] = (-1*S[i]*pow((((double)XRANGE/Px)*((double)YRANGE/Py)),2)+(T[i-1]+T[i+1])*pow(((double)YRANGE/Py),2)+ \
					(T[i-Px]+T[i+Px])*pow(((double)XRANGE/Px),2))/(2*pow(((double)XRANGE/Px),2)+2*pow(((double)YRANGE/Py),2));
					if(unlock == 1) {
						omp_unset_lock(&mutexes[mutex_map(i, id)]);
						unlock = 0;
					}

				}
			}
	}

	double Xmin = 0;
	double Ymin = 0;
	double Xmax = XRANGE;
	double Ymax = YRANGE;
	VTK_out(Px, Py, &Xmin, &Xmax, &Ymin, &Ymax, T, 0);
	print_grid(T, grid_size, Px);




	free(S);
	free(T);
	return 0;
}
