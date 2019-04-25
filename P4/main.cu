#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define XRANGE 1
#define YRANGE 1

#define THRESH 1e-12

extern "C"
void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);

//__device__ volatile int d_error_mutex = 0;
//__device__ double * d_conv_error;

__global__
void prepare_grids(double *T, double * T_tmp, double *S, double * errors, long * grid_size, long * internal_size, int * Px, int * Py) {
	long id = blockIdx.x*blockDim.x + threadIdx.x;
	long mapped_id = id-(2*(id/Px[0])-1)-(Px[0]);
	if(id >= grid_size[0])
		return;
	double val = (id%Px[0])*((double)XRANGE/Px[0])*powf(2.718281828, (id/Px[0])*((double)YRANGE/Py[0]));
	S[id] = val;
	if(id/Px[0]==0 || id/Px[0]==Py[0]-1 || id%Px[0]==0 || id%Px[0]==Px[0]-1) { 
		T[id] = val;
		T_tmp[id] = val;
		return;
	} else {
		T[id] = 0;
		T_tmp[id] = 0;
	}
	errors[mapped_id] = 0;
}

__global__
void update_temporary(double * T, double * T_tmp, double * S, double * errors, long * grid_size, int * Px, int * Py) {
	long id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id>=grid_size[0] || id/Px[0]==0 || id/Px[0]==Py[0]-1 || id%Px[0]==0 || id%Px[0]==Px[0]-1)
		return;
	T_tmp[id] = (-1*S[id]*powf((((double)XRANGE/Px[0])*((double)YRANGE/Py[0])),2)+(T[id-1]+T[id+1])*powf(((double)YRANGE/Py[0]),2)+ \
	(T[id-Px[0]]+T[id+Px[0]])*powf(((double)XRANGE/Px[0]),2))/(2*powf(((double)XRANGE/Px[0]),2)+2*powf(((double)YRANGE/Py[0]),2));
}

__global__
void update_real(double * T, double * T_tmp, double * S, double * errors, long * grid_size, int * Px, int * Py) {
	long id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id>=grid_size[0] || id/Px[0]==0 || id/Px[0]==Py[0]-1 || id%Px[0]==0 || id%Px[0]==Px[0]-1)
		return;
	T[id] = T_tmp[id];
}

__global__ 
void get_abs_error(double *T, double *S, long * grid_size, double * d_abs_error) {	
	long id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id >= grid_size[0])
		return;
	if(id==0)
		*d_abs_error = 0;
	__syncthreads();
	double val = fabs(T[id] - S[id]);
	val = powf(val,3);
	atomicAdd(d_abs_error, val);
}

__global__
void get_error(double *T, double *T_tmp, long * grid_size, int * Px, int * Py, double * d_conv_error) {
	long id = blockIdx.x*blockDim.x + threadIdx.x;
	if(id >= grid_size[0])
		return;
	if(id==0)
		*d_conv_error = 0;
	__syncthreads();
	double val = fabs(T[id] - T_tmp[id]);
	val = powf(val,3);
	atomicAdd(d_conv_error, val);
}

int main(int argc, char **argv) {
	if(argc!=3) {
		printf("Usage:\n\t./main [xdim] [ydim]\n");
		exit(0);
	}
	int h_Px = atoi(argv[1]);
	int h_Py = atoi(argv[2]);
	long h_grid_size = h_Px*h_Py;
	long h_internal_size = h_grid_size - 2*h_Px - 2*h_Py + 4;
	int *d_Px, *d_Py;
	long *d_grid_size, *d_internal_size;
	cudaMalloc((void**)&d_Px, sizeof(int));
	cudaMalloc((void**)&d_Py, sizeof(int));
	cudaMalloc((void**)&d_grid_size, sizeof(long));
	cudaMalloc((void**)&d_internal_size, sizeof(long));
	cudaMemcpy(d_Px, &h_Px, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Py, &h_Py, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_grid_size, &h_grid_size, sizeof(long), cudaMemcpyHostToDevice);
	cudaMemcpy(d_internal_size, &h_internal_size, sizeof(long), cudaMemcpyHostToDevice);


	// Allocate grids in Host to recieve from GPU (pinned memory)
	double * h_S;
	double * h_T;
	double * h_T_tmp;
	cudaMallocHost((void**)&h_S, h_grid_size*sizeof(double));
	cudaMallocHost((void**)&h_T, h_grid_size*sizeof(double));	
	cudaMallocHost((void**)&h_T_tmp, h_grid_size*sizeof(double));

	// Allocate grids in GPU memory
	double * d_S;
	double * d_T;
	double * d_T_tmp;
	double * d_errors;
	cudaMalloc((double**)&d_S, h_grid_size*sizeof(double));
	cudaMalloc((double**)&d_T, h_grid_size*sizeof(double));	
	cudaMalloc((double**)&d_T_tmp, h_grid_size*sizeof(double));
	cudaMalloc((double**)&d_errors, h_internal_size*sizeof(double));
	double * d_conv_error;
	double * d_abs_error;
	cudaMalloc((double**)&d_conv_error, sizeof(double));
	cudaMalloc((double**)&d_abs_error, sizeof(double));

	int blocks = (ceil((double)(h_Px*h_Py)/1000));
	int threadsperblock = 1000;
	printf("Running on %d blocks each with %d threads\n",blocks,threadsperblock);

	prepare_grids<<<blocks, threadsperblock>>>(d_T,d_T_tmp,d_S,d_errors,d_grid_size,d_internal_size,d_Px,d_Py);
	cudaDeviceSynchronize();
	int iter = 0;
	double h_conv_error = THRESH+1;
	//double h_abs_error = 0;
	while(h_conv_error > THRESH) {
		update_temporary<<<blocks,threadsperblock>>>(d_T,d_T_tmp,d_S,d_errors,d_grid_size,d_Px,d_Py);
		cudaDeviceSynchronize();
		if(iter%1000==0) {
			//get_error<<<blocks,threadsperblock>>>(d_T,d_T_tmp,d_grid_size,d_Px,d_Py,d_conv_error);
			//get_abs_error<<<blocks,threadsperblock>>>(d_T,d_S,d_grid_size,d_abs_error);
			//cudaDeviceSynchronize();
			//cudaMemcpy(&h_conv_error, d_conv_error, sizeof(double), cudaMemcpyDeviceToHost);
			//cudaMemcpy(&h_abs_error, d_abs_error, sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_T, d_T, h_grid_size*sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(h_T_tmp, d_T_tmp, h_grid_size*sizeof(double), cudaMemcpyDeviceToHost);
			h_conv_error = 0;
			for(long i=0; i<h_grid_size; i++) {
				double val = fabs(h_T_tmp[i] - h_T[i]);
				if(val > h_conv_error)
					h_conv_error = val;
			}		
			//h_abs_error = powf(h_abs_error, 0.333);
			//h_conv_error = powf(h_conv_error, 0.333);
			printf("iter = %d... conv Error = %.10e\n", iter, h_conv_error);
		}
		update_real<<<blocks,threadsperblock>>>(d_T,d_T_tmp,d_S,d_errors,d_grid_size,d_Px,d_Py);
		iter++;
	}
	printf("Finished\n");
	
	cudaMemcpy(h_T, d_T, h_grid_size*sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_S, d_S, h_grid_size*sizeof(double), cudaMemcpyDeviceToHost);
	double max = 0;
	for(long i=0; i<h_grid_size; i++) {
		double val = fabs(h_S[i] - h_T[i]);
		if(val > max)
			max = val;
	}
	printf("ABSLOUTE ERROR = %.10e\n", max);
	// Output .vtk file for ParaView
	double Xmin = 0;
	double Ymin = 0;
	double Xmax = XRANGE;
	double Ymax = YRANGE;
	VTK_out(h_Px, h_Py, &Xmin, &Xmax, &Ymin, &Ymax, h_T, 0);
	
}
