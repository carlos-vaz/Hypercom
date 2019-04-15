#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define XRANGE 2
#define YRANGE 1

#define TOLERANCE	1e-12
#define ROUNDS 		1000

int Px, Py, Tx, Ty, Ptx, Pty, shared_pts;
long grid_size;
double conv_error;


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

double AbsoluteError(double * T, double * S) {
	double max = 0, this = 0;
	for(int i=0; i<grid_size; i++) {
		this = fabs(T[i]-S[i]);
		if(max > this) {
			max = this;
		}
	}
	return max;
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
	int per_thread = Ptx*Pty;
	shared_pts = (Tx-1)*(Py-2) + (Ty-1)*(Px-2) - 4*(Tx-1)*(Ty-1);
	if((Px%Tx!=0) || (Py%Ty!=0)) {
		printf("Threads in X dimension must divide Points in X dimension, and likewise for Y dimension\n");
		exit(-2);
	}
	printf("Running with parameters:\n\tPoints X = %d, Threads X = %d\n\tPoints Y = %d, Threads Y = %d\n", Px, Tx, Py, Ty);

	// Make Grids
	double * S = (double * )malloc(grid_size*sizeof(double));
	double * T = (double * )calloc(grid_size, sizeof(double));
	double * T_tmp = (double * )calloc(grid_size, sizeof(double));
	for(long i=0; i<grid_size; i++) {
		double val = (i%Px)*((double)XRANGE/Px)*powf(2.718281828, (i/Px)*((double)YRANGE/Py));
		//printf("TEST: %ld/%d = %lf\n", i, Px, (double)(i/Px));
		S[i] = val;
		if(i/Px==0 || i/Px==Py-1 || i%Px==0 || i%Px==Px-1) {
			T[i] = val;
			T_tmp[i] = val;
		}
	}

	//print_grid(T, grid_size, Px);

	omp_set_num_threads(Tx*Ty);
	int id, i, unlock=0;
	int count = 0;
	double myAbsError = 0, AbsError = 0;
	double conv_error = TOLERANCE+1, mymax, candidate;
	while(conv_error > TOLERANCE) {

		#pragma omp parallel default(none) shared(S,T,T_tmp,grid_size,Tx,Ty,Px,Py,Ptx,Pty,conv_error,count,per_thread,AbsError) private(id,i,unlock,mymax,myAbsError,candidate)
		{
			id = omp_get_thread_num();
			//printf("Hello from thread %d\n", id);

			#pragma omp for
				for(i=0; i<grid_size; i++) {
					//printf("\tSubroutine from thread %d (i=%d, Px = %d)\n", id, i, Px);
					if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1) {
						T_tmp[i] = (-1*S[i]*pow((((double)XRANGE/Px)*((double)YRANGE/Py)),2)+(T[i-1]+T[i+1])*pow(((double)YRANGE/Py),2)+ \
						(T[i-Px]+T[i+Px])*pow(((double)XRANGE/Px),2))/(2*pow(((double)XRANGE/Px),2)+2*pow(((double)YRANGE/Py),2));
					}
				}

			if(count%ROUNDS==0) {
				mymax = 0;
				conv_error = 0;
				#pragma omp for
				for(i=0; i<grid_size; i++) {
					if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1) {
						candidate = fabs(T_tmp[i] - T[i]);
						//printf("Candidate: %lf\n", candidate);
						if(mymax < candidate)
							mymax = candidate;
					}
				}
				#pragma omp critical
				{
					if(mymax > conv_error)
						conv_error = mymax;
				}


				#pragma omp for
				for(i=0; i<grid_size; i++) {
					if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1) {
						candidate = fabs(S[i] - T[i]);
						//printf("Candidate: %lf\n", candidate);
						if(myAbsError < candidate)
							myAbsError = candidate;
					}
				}
				#pragma omp critical
				{
					if(myAbsError > AbsError)
						AbsError = myAbsError;
				}
				
		
			}

			//printf("(thread %d): mymax = %lf, conv_error = %lf\n", mymax, conv_error);


			//printf("\t (thread %d AFTER):  CONV ERROR= %lf\n", conv_error);

			#pragma omp barrier
			if(id==0 & count%ROUNDS==0) {
				printf("%d: CONV Error this round = %.10e\n", count, conv_error);
				printf("%d: ABS Error this round = %.10e\n", count, AbsError);
			}


			/*#pragma omp for
				for(i=0; i<grid_size; i++) {
					if(i/Px!=0 && i/Px!=Py-1 && i%Px!=0 && i%Px!=Px-1)
						T[i] = T_tmp[i];
				}*/
			memcpy(&T[id*per_thread], &T_tmp[id*per_thread], per_thread*sizeof(double));

			
		}
		//printf("Left parallel part\n");
		count++;
	}


	double Xmin = 0;
	double Ymin = 0;
	double Xmax = XRANGE;
	double Ymax = YRANGE;
	VTK_out(Px, Py, &Xmin, &Xmax, &Ymin, &Ymax, T, 0);
	//print_grid(T, grid_size, Px);




	free(S);
	free(T);
	free(T_tmp);
	return 0;
}
