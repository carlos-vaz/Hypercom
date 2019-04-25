#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define N 1000000
#define THREADS 1024

__global__
void kernel(double * F, double * final, int threads) {
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	int chunksize = ceil((double)N/THREADS);
	for(int i=id*chunksize; i<min(id*chunksize+chunksize, N); i++) {
		F[i] = (double)4/(1+powf(((double)i/N),2));
	}
	__syncthreads();
	for(int i=id*chunksize; i<min(id*chunksize+chunksize, N); i++) {
		atomicAdd(final, (F[i] + F[i+1])/(2*N));
	}
}

int main(int argc, char **argv) {
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);	
	double * d_F;
	double * d_value;
	cudaMalloc((void**)&d_F, N*sizeof(double));	
	cudaMalloc((void**)&d_value, sizeof(double));
	double h_value;

	cudaEventRecord(start);
	kernel<<<1,THREADS>>>(d_F,d_value,THREADS);
	cudaEventRecord(stop);

	cudaMemcpy(&h_value, d_value, sizeof(double), cudaMemcpyDeviceToHost);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("VALUE:\t%lf\n", h_value);
	printf("ELAPSED:\t%f\n", milliseconds);
}
