#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define XRANGE 1
#define YRANGE 1

#define THRESH 1e-12

__device__
double reduce_conv_error(double * T, double * T_tmp, double * errors, int internal_size, int map_id) {
	int stride = 2;
	while(stride < internal_size) {
		if(map_id % stride == 0) {
			errors[map_id] = fmax(errors[map_id], errors[map_id+stride/2]);
			if(internal_size % stride == (stride/2) && map_id==internal_size-3*stride/2) {
				// This threads should collect the max of THREE values, not two
				errors[map_id] = fmax(errors[map_id], errors[map_id+stride]);
			}
		}
		stride <<= 1;
		__syncthreads();
	}
	return errors[0];
}

__global__
void kernel(double * T, double * T_tmp, double * S, double * errors, int grid_size, int internal_size, int Px, int Py) {
	// First, fill grids with boundary conditions
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id >= grid_size)
		return;
	double val = (id%Px)*((double)XRANGE/Px)*powf(2.718281828, (id/Px)*((double)YRANGE/Py));
	S[id] = val;
	T[id] = 0;
	T_tmp[id] = 0;
	if(id/Px==0 || id/Px==Py-1 || id%Px==0 || id%Px==Px-1) {
		T[id] = S[id];
		T_tmp[id] = S[id];
		return;
	}

	__syncthreads();

	// Then, begin computing solution
	int iter = 0;
	double error = THRESH+1;
	int map_id;
	while(iter < 1000) {
		T_tmp[id] = (-1*S[id]*pow((((double)XRANGE/Px)*((double)YRANGE/Py)),2)+(T[id-1]+T[id+1])*pow(((double)YRANGE/Py),2)+ \
		(T[id-Px]+T[id+Px])*pow(((double)XRANGE/Px),2))/(2*pow(((double)XRANGE/Px),2)+2*pow(((double)YRANGE/Py),2));
		__syncthreads();
		T[id] = T_tmp[id];
		__syncthreads();
		if(iter%1000==0) {
			// Remap id to index only internal grid points
			map_id = id-(2*(id/Px)-1)-(Px);
			error = reduce_conv_error(T, T_tmp, errors, internal_size, map_id);
		}
		iter++;
	}
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
	long internal_size = grid_size - 2*Px - 2*Py + 4;

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
	double * d_errors;
	cudaMalloc((double**)&d_S, grid_size*sizeof(double));
	cudaMalloc((double**)&d_T, grid_size*sizeof(double));	
	cudaMalloc((double**)&d_T_tmp, grid_size*sizeof(double));
	cudaMalloc((double**)&d_errors, internal_size*sizeof(double));

	int blocks = ceil((double)Px*Py/1000);
	int threadsperblock = 1000;
	printf("Running on %d blocks each with %d threads\n",blocks,threadsperblock);
	kernel<<<blocks,threadsperblock>>>(d_T,d_T_tmp,d_S,d_errors,grid_size,internal_size,Px,Py);

	cudaDeviceSynchronize();

	cudaMemcpy(h_T, d_T, grid_size*sizeof(double), cudaMemcpyDeviceToHost);
	for(int i=0; i<grid_size; i++) {
		printf("%lf ", h_T[i]);
		if(i%Px==Px-1)
			printf("...\n");
	}

}
