#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define XRANGE 1
#define YRANGE 1

__global__
void kernel(double * T, double * T_tmp, double * S, int grid_size, int Px, int Py) {
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id >= grid_size)
		return;
	double val = (id%Px)*((double)XRANGE/Px)*powf(2.718281828, (id/Px)*((double)YRANGE/Py));
	S[id] = val;
	T[id] = 0;
	T_tmp[id] = 0;
	if(id/Px!=0 && id/Px!=Py-1 && id%Px!=0 && id%Px!=Px-1)
		return;
	T[id] 	  = S[id];	
	T_tmp[id] = S[id];
	__syncthreads();


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

	// Allocate grids in Host to recieve from GPU (pinned memory)
	double * h_S;
	double * h_T;
	double * h_T_tmp;
	cudaMallocHost((void**)&h_S, grid_size*sizeof(double));
	cudaMallocHost((void**)&h_T, grid_size*sizeof(double));
	cudaMallocHost((void**)&h_T_tmp, grid_size*sizeof(double));
	
	// Allocate grids in GPU memory
	double * d_S;
	double * d_T;
	double * d_T_tmp;
	cudaMalloc((double**)&d_S, grid_size*sizeof(double));
	cudaMalloc((double**)&d_T, grid_size*sizeof(double));	
	cudaMalloc((double**)&d_T_tmp, grid_size*sizeof(double));
	//cudaMemcpy(d_S, h_S, grid_size*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_T, h_T_tmp, grid_size*sizeof(double), cudaMemcpyHostToDevice);
	//cudaMemcpy(d_T_tmp, h_T_tmp, grid_size*sizeof(double), cudaMemcpyHostToDevice);

	int blocks = ceil((double)Px*Py/1000);
	int threadsperblock = 1000;
	printf("Running on %d blocks each with %d threads\n",blocks,threadsperblock);
	kernel<<<blocks,threadsperblock>>>(d_T,d_T_tmp,d_S,grid_size,Px,Py);

	cudaDeviceSynchronize();

	cudaMemcpy(h_T, d_T, grid_size*sizeof(double), cudaMemcpyDeviceToHost);
	for(int i=0; i<grid_size; i++) {
		printf("%lf ", h_T[i]);
		if(i%Px==Px-1)
			printf("...\n");
	}

}
