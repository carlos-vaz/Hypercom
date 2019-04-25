#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define N 1000000
#define THREADS 512

__global__
void kernel(double * F, int threads) {
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id >= N)
		return;
	if(id%THREADS==0)
		for(int i=0; i<THREADS; i++) {
			F[id + i] = (double)4/(1+powf(((double)id/N),2));
		}
	__syncthreads();
	
	
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
}
