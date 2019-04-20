#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define XRANGE 1
#define YRANGE 1

__global__
void kernel(double * T, double * S) {
	int id = blockIdx.x*blockDim.x + threadIdx.x;
}

int Px, Py, grid_size;

int main(int argc, char **argv) {
	if(argc!=3) {
		printf("Usage:\n\t./main [xdim] [ydim]\n");
		exit(0);
	}
	Px = atoi(argv[1]);
	Py = atoi(argv[2]);
	long grid_size = Px*Py;

	// Make initial grids (pinned memory)
	double * h_S;
	double * h_T;
	double * h_T_tmp;
	cudaMallocHost((void**)&h_S, grid_size*sizeof(double));
	cudaMallocHost((void**)&h_T, grid_size*sizeof(double));
	cudaMallocHost((void**)&h_T_tmp, grid_size*sizeof(double));
	for(long i=0; i<grid_size; i++) {
		double val = (i%Px)*((double)XRANGE/Px)*powf(2.718281828, (i/Px)*((double)YRANGE/Py));
		h_S[i] = val;
		if(i/Px==0 || i/Px==Py-1 || i%Px==0 || i%Px==Px-1) {
			h_T[i] = val;
			h_T_tmp[i] = val;
		}
	}
	
	// Transfer grids to GPU memory
	double * d_S;
	double * d_T;
	double * d_T_tmp;
	cudaMalloc((double**)&d_S, grid_size*sizeof(double));
	cudaMalloc((double**)&d_T, grid_size*sizeof(double));	
	cudaMalloc((double**)&d_T_tmp, grid_size*sizeof(double));
	cudaMemcpy(d_S, h_S, grid_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_T, h_T_tmp, grid_size*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_T_tmp, h_T_tmp, grid_size*sizeof(double), cudaMemcpyHostToDevice);

	int blocks = 1;
	int threadsperblock = 1;
	kernel<<<blocks,threadsperblock>>>(d_T,d_S);

}
