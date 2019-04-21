#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"

#define XRANGE 1
#define YRANGE 1

#define THRESH 1e-12

extern 
void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);

__device__
double reduce_conv_error(double * T, double * T_tmp, double * errors, long internal_size, int id, int map_id) {
	errors[map_id] = fabs(T[id] - T_tmp[id]);

	int stride = 2;
	int num_active = internal_size/2;
	int prev_active = internal_size;
	int  active = map_id%2==0 ? 1 : 0;
	if(internal_size % 2 == 1 && map_id == internal_size-1) active = 0;
	while(num_active > 0) {
		if(active==1) {
			errors[map_id] = fmax(errors[map_id], errors[map_id + stride/2]);
			//printf("[%d] <-- [%d]\n", map_id, map_id+stride/2);
			if(prev_active%2==1 && map_id==(num_active-1)*stride) {
				errors[map_id] = fmax(errors[map_id], errors[map_id + stride]);
				//printf("special. [%d] <-- [%d]\n", map_id, map_id+stride);//, map_id, map_id + stride);
			}		
		}
		stride <<= 1;
		if(map_id % stride != 0 || (num_active%2==1 && map_id==(num_active-1)*stride/2))
			active = 0;
		prev_active = num_active;
		num_active /= 2;
		//if(map_id==0) printf("----------\n");
		__syncthreads();
	}
	return errors[0];
}


__global__
void kernel(double * T, double * T_tmp, double * S, double * errors, long grid_size, long internal_size, int Px, int Py) {
	// First, fill grids with boundary conditions
	int id = blockIdx.x*blockDim.x + threadIdx.x;
	//int active = 1;
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
	int mapped_id = id-(2*(id/Px)-1)-(Px);
	int iter = 0;
	double error = THRESH+1;
	while(error > THRESH) {
		T_tmp[id] = (-1*S[id]*pow((((double)XRANGE/Px)*((double)YRANGE/Py)),2)+(T[id-1]+T[id+1])*pow(((double)YRANGE/Py),2)+ \
		(T[id-Px]+T[id+Px])*pow(((double)XRANGE/Px),2))/(2*pow(((double)XRANGE/Px),2)+2*pow(((double)YRANGE/Py),2));
		__syncthreads();
		if(iter%1000==0) {
			//errors[mapped_id] = 1; // test reduce to sum
			error = reduce_conv_error(T, T_tmp, errors, internal_size, id, mapped_id);
			if(mapped_id==0) printf("iter %d Error = %.10e\n", iter, error);
		}
		T[id] = T_tmp[id];
		__syncthreads();
		iter++;
	}
}

int Px, Py;

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

	// Output .vtk file for ParaView
	double Xmin = 0;
	double Ymin = 0;
	double Xmax = XRANGE;
	double Ymax = YRANGE;
	VTK_out(Px, Py, &Xmin, &Xmax, &Ymin, &Ymax, T, 0);

}
