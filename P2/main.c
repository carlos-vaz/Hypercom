#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

#define ERROR_THRESH  	1e-4	// Max tolerable error (infinity norm)
#define ANNOUNCER_PROC 	0
//#define	   index(x,y)	(x+y*proc_pts[0])
#define	   up(i)	(i+proc_pts[0])
#define	   down(i)	(i-proc_pts[0])
#define    right(i)	(i+1)
#define    left(i)	(i-1)
	

extern 
void VTK_out(const int N, const int M, const double *Xmin, const double *Xmax,
             const double *Ymin, const double *Ymax, const double *T,
             const int index);

/*
 * Return value:
	0 = Fail
	1 = Pass
 */
int verify_error(double *T, double *T_prev, int len) {
	for(int i=0; i<len; i++) {
		if(fabs(T[i]-T_prev[i])>ERROR_THRESH) {
			return 0;
		}
	}
	return 1;
}

double get_error(double *T, double *T_prev, int len) {
	double ret = 0;
	double tmp;
	for(int i=0; i<len; i++) {
		tmp = fabs(T[i]-T_prev[i]);
		if(tmp>ret) {
			ret = tmp;
		}
	}
	printf("\tget_error returnint %lf\n", ret);
	return ret;
}

int myrank, rank_2d, mycoord[2], np, dims_procs[2], num_points, dims_pts[2], proc_pts[2], proc_size, \
	ranks_around[4] = {-1,-1,-1,-1}; // {right, left, up, down}
double ranges[2], deltas[2];


/*
 * Each process is assigned a grid section via MPI_Cart_create, and the corresponding 
 * region of the data file is read into each process using MPI_Type_vector and 
 * MPI_File_read_all. 
 *  
 */
int main(int argc, char* argv[]) {
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	/*
	 * Check for correct number of arguments
	 * Arguments in the form: Y procs, X procs (MATRIX FORM)
	 */
	if(argc!=4) {
		if(myrank == ANNOUNCER_PROC)
			printf("ERROR: arguments to main: mpirun -np [#procs] main [sourcefile] [# procs Y] [# procs X] \n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}


	/*
	 * Read the file metadata and check for do-ablility
	 * File metadata in the form: X dimensions, Y dimensions, X range, Y range (CARTESIAN FORM)
	 */
	dims_procs[0] = atoi(argv[3]);
	dims_procs[1] = atoi(argv[2]);
	if(myrank == ANNOUNCER_PROC) printf("Will solve %s with %d (x) by %d (y) processes\n", argv[1], dims_procs[0], dims_procs[1]);
	MPI_File file; 
	MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
	MPI_File_read_all(file, &dims_pts, 2, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read_at_all(file, 2*sizeof(int), &ranges, 2, MPI_DOUBLE, MPI_STATUS_IGNORE);
	deltas[0] = ranges[0]/(double)dims_pts[0];
	deltas[1] = ranges[1]/(double)dims_pts[1];
	if(myrank == ANNOUNCER_PROC) printf("FILE READ... DIMENSIONS (MATRIX FORM): %d BY %d\n", dims_pts[1], dims_pts[0]);
	if( (dims_pts[0]%dims_procs[0]!=0 || dims_pts[1]%dims_procs[1]!=0) && myrank==ANNOUNCER_PROC) {
		printf("Number of points must be divisible by number of procs in that dimension.\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}
	if(np!=dims_procs[0]*dims_procs[1] && myrank == ANNOUNCER_PROC) {
		printf("Allocated %d processes, but user specified %d*%d = %d processes.\n", \
				np, dims_procs[1], dims_procs[0], dims_procs[1]*dims_procs[0]);
		MPI_Abort(MPI_COMM_WORLD, -1);	
	}
	MPI_Barrier(MPI_COMM_WORLD); // wait in case proc 0 aborts


	/*
	 * Create a cartesian topology. This lets us identify processes by their 
	 * coordinates instead of their ranks. 
	 */
	MPI_Comm comm2d;
	int periodic[2] = {0,0};	// What is periodic?
	MPI_Cart_create(MPI_COMM_WORLD, 2,dims_procs,periodic,1,&comm2d);
	MPI_Cart_coords(comm2d,myrank, 2,mycoord);
	MPI_Cart_rank(comm2d,mycoord,&rank_2d);
	printf("I am %d: (%d,%d); originally %d\n",rank_2d,mycoord[0],mycoord[1],myrank);
	


	/*
	 * Read the file contents into vector
	 */
	MPI_Datatype vector; 
	proc_pts[0] = dims_pts[0]/dims_procs[0]; // # pts in Y dim of process partition
	proc_pts[1] = dims_pts[1]/dims_procs[1]; // # pts in X dim of process partition
	proc_size = proc_pts[0]*proc_pts[1];
	double *v = (double*)malloc(proc_size*sizeof(double));
	if( myrank == ANNOUNCER_PROC) printf("# blocks per process = %d\t# pts per proc (X) = %d\t# pts per proc (Y) = %d\n", \
							proc_pts[1], proc_pts[0], proc_pts[1]);
	MPI_Type_vector(proc_pts[1], proc_pts[0], dims_pts[0], MPI_DOUBLE, &vector);
	MPI_Type_commit(&vector);	
	MPI_File_set_view(file, 2*sizeof(int)+2*sizeof(double)+(mycoord[0]*proc_pts[0]+mycoord[1]*proc_size*dims_procs[0])*sizeof(double), \
								MPI_DOUBLE, vector, "native", MPI_INFO_NULL);	
	MPI_File_read_all(file, v, proc_size, MPI_DOUBLE, MPI_STATUS_IGNORE);

/*	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0; i<np; i++) {
		if(myrank==i) {
			printf("\nRANK %d (%d, %d)\n", myrank, mycoord[0], mycoord[1]);
			for(int j=0; j<proc_size; j++) {
				printf("%d, ", (int)v[j]);
			}
			printf("\n\n");
		}
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
	}
*/
	/*
	 * Each process creates a temperature vector T and fills it with 
	 * 0, except for at the boundaries, which take the initial value 
	 * of xe^y (so we copy from vector 'v' at boundary points)
	 */
	int bound_south = 0;
	int bound_north = 0;
	int bound_east = 0;
	int bound_west = 0;
	double *T = (double*)calloc(proc_size, sizeof(double));
	if(mycoord[1]==0) {
		//printf("Process (%d, %d): I have a southern boundary\n", mycoord[0], mycoord[1]);
		bound_south = 1;
		memcpy(v, T, proc_pts[0]*sizeof(double));
	}
	if(mycoord[1]==dims_procs[1]-1) {
		//printf("Process (%d, %d): I have a northern boundary\n", mycoord[0], mycoord[1]);
		bound_north = 1;
		for(int x=0; x<proc_pts[0]; x++)
			T[(proc_pts[1]-1)*proc_pts[0]+x] = v[(proc_pts[1]-1)*proc_pts[0]+x]; // copy from 'v'
	}
	if(mycoord[0]==0) {
		//printf("Process (%d, %d): I have a western boundary\n", mycoord[0], mycoord[1]);
		bound_west = 1;
		for(int y=0; y<proc_pts[1]; y++)
			T[y*proc_pts[0]] = v[y*proc_pts[0]]; // copy from 'v'
		
	}
	if(mycoord[0]==dims_procs[0]-1) {
		//printf("Process (%d, %d): I have an eastern boundary\n", mycoord[0], mycoord[1]);
		bound_east = 1;
		for(int y=0; y<proc_pts[1]; y++)
			T[(y+1)*proc_pts[0]-1] = v[(y+1)*proc_pts[0]-1];  // copy from 'v'

	}


	// Preparation complete... Start the timer

	// Allocate 1 extra element to signal when the neighbor should remember
	// the sent buffer because the sender will stop sending (due to exit)
	double *send_south = (double*)calloc((proc_pts[0]+1), sizeof(double));
	double *send_north = (double*)calloc((proc_pts[0]+1), sizeof(double));
	double *send_east  = (double*)calloc((proc_pts[1]+1), sizeof(double));
	double *send_west  = (double*)calloc((proc_pts[1]+1), sizeof(double));
	double *recv_south = (double*)calloc((proc_pts[0]+1), sizeof(double));
	double *recv_north = (double*)calloc((proc_pts[0]+1), sizeof(double));
	double *recv_east  = (double*)calloc((proc_pts[1]+1), sizeof(double));
	double *recv_west  = (double*)calloc((proc_pts[1]+1), sizeof(double));
	int got_south = 0;
	int got_north = 0;
	int got_east  = 0;
	int got_west  = 0;
	MPI_Cart_shift(comm2d, 0, +1, &rank_2d, &ranks_around[0]);
	MPI_Cart_shift(comm2d, 0, -1, &rank_2d, &ranks_around[1]);
	MPI_Cart_shift(comm2d, 1, +1, &rank_2d, &ranks_around[2]);
	MPI_Cart_shift(comm2d, 1, -1, &rank_2d, &ranks_around[3]);

	int remember_south = 0;
	int remember_north = 0;
	int remember_east  = 0;
	int remember_west  = 0;

	// Rank order: east/west/north/south = 0/1/2/3
	MPI_Request req[8];
	MPI_Status stati[8];
	double * test_buffer = (double*)malloc(proc_size*sizeof(double));
	int count = 0;
	int will_break = 0;
	
	while(count < 10000) {
//	while(1) {
		/*
		 * Post a non-blocking send and a non-blocking receive to all neighbors.
		 * While you update your internal temperatures, hopefully the requests
		 * will go through. 
		 */
		if(bound_south==0 && remember_south==0) {
			got_south = 1;
			MPI_Irecv(recv_south, proc_pts[0]+1, MPI_DOUBLE, ranks_around[3] /*southern rank*/ \
									, 2 /*northernly tag*/, comm2d, &req[0]);
			if(recv_south[proc_pts[0]]==1) // did Southern neighbr signal to me to remember recv_south?
				remember_south = 1;
			//memcpy(send_south, T, proc_pts[0]*sizeof(double)); // Why copy if T[0->xdim] won't change?
			for(int i=0; i<proc_pts[0]; i++)
				send_south[i] = T[i];
			if(will_break)
				send_south[proc_pts[0]] = 1; // signal to Southern neighbor to remember this buffer
			MPI_Isend(send_south, proc_pts[0]+1, MPI_DOUBLE, ranks_around[3] /*southern rank*/, 3/*southernly tag*/, comm2d, &req[1]);
		}
		if(bound_north==0 && remember_north==0) {
			got_north = 1;
			//memcpy(send_north, &T[proc_size-proc_pts[0]], proc_pts[0]*sizeof(double)); // Why copy if T[0->xdim] won't change?
			for(int i=0; i<proc_pts[0]; i++)
				send_north[i] = T[proc_size-proc_pts[0]+i];
			if(will_break)
				send_north[proc_pts[0]] = 1; // signal to Northern neighbor to remember this buffer

			MPI_Isend(send_north, proc_pts[0]+1, MPI_DOUBLE, ranks_around[2] /*northern rank*/, 2/*northernly tag*/, comm2d, &req[2]);
			MPI_Irecv(recv_north, proc_pts[0]+1, MPI_DOUBLE, ranks_around[2] /*northern rank*/ \
									, 3 /*southernly tag*/, comm2d, &req[3]);
			if(recv_north[proc_pts[0]]==1) // did Northern neighbr signal to me to remember recv_north?
				remember_north = 1;
		}
		if(bound_east==0 && remember_east==0) {
			got_east = 1;
			MPI_Irecv(recv_east, proc_pts[1]+1, MPI_DOUBLE, ranks_around[0] /*eastern rank*/ \
									, 1 /*westernly tag*/, comm2d, &req[4]);
			if(recv_east[proc_pts[1]]==1) // did Eastern neighbr signal to me to remember recv_east?
				remember_east = 1;
			//got_east=1;
			// Copy eastern buffer to send_east
			for(int i=0; i<proc_pts[1]; i++)
				send_east[i] = T[proc_pts[0]-1+i*proc_pts[0]];
			if(will_break)
				send_east[proc_pts[1]] = 1; // signal to Eastern neighbor to remember this buffer
			MPI_Isend(send_east, proc_pts[1]+1, MPI_DOUBLE, ranks_around[0] /*eastern rank*/, 0/*easternly tag*/, comm2d, &req[5]);
		}
		if(bound_west==0 && remember_west==0) {
			got_west = 1;
			// Copy western buffer to send_west
			for(int i=0; i<proc_pts[1]; i++)
				send_west[i] = T[i*proc_pts[0]];
			if(will_break)
				send_west[proc_pts[1]] = 1; // signal to Western neighbor to remember this buffer
			MPI_Isend(send_west, proc_pts[1]+1, MPI_DOUBLE, ranks_around[1] /*western rank*/, 1/*westernly tag*/, comm2d, &req[6]);
			MPI_Irecv(recv_west, proc_pts[1]+1, MPI_DOUBLE, ranks_around[1] /*western rank*/ \
									, 0 /*easternly tag*/, comm2d, &req[7]);
			if(recv_west[proc_pts[1]]==1) // did Western neighbr signal to me to remember recv_west?
				remember_west = 1;
			//got_west=1;
		}


		/*
		 * Compute & update temperatures on block interiors (Gauss-Seidel + block Jacobi) for
		 * 1 iteration. 
		 */
		int i=0, x=0, y=0;
		if(got_south==1 && got_west==1) {
			//i = index(0,0);
			i = 0;
			T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(recv_west[0]+T[right(i)])*pow(deltas[1],2)+ \
				(recv_south[0]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
		}
		if(got_south==1) {
			for(x=1; x<proc_pts[0]-1; x++) {
				//i = index(x,0);
				i = x;
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+T[right(i)])*pow(deltas[1],2)+ \
					(recv_south[x]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
		}
		if(got_south==1 && got_east==1) {
			//i = index(proc_pts[0]-1,0);
			i = proc_pts[0]-1;
			T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+recv_east[0])*pow(deltas[1],2)+ \
				(recv_south[proc_pts[0]-1]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
		}
		for(y=1; y<proc_pts[1]-1; y++) {
			if(got_west==1) {
				//i = index(0,y);
				i = y*proc_pts[0];
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(recv_west[y]+T[right(i)])*pow(deltas[1],2)+ \
					(T[down(i)]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
			for(x=1; x<proc_pts[0]-1; x++) {
				//i = index(x,y);
				i = x+y*proc_pts[0];
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+T[right(i)])*pow(deltas[1],2)+ \
					(T[down(i)]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
			if(got_east==1) {
				//i = index(proc_pts[0]-1, y);
				i = proc_pts[0]-1+y*proc_pts[0];
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+recv_east[y])*pow(deltas[1],2)+ \
					(T[down(i)]+T[up(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
		}
		if(got_north==1 && got_west==1) {
			//i = index(0,proc_pts[1]-1);
			i = (proc_pts[1]-1)*proc_pts[0];
			T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(recv_west[proc_pts[1]-1]+T[right(i)])*pow(deltas[1],2)+ \
				(recv_north[0]+T[down(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
		}
		if(got_north==1) {
			for(x=1; x<proc_pts[0]-1; x++) {
				//i = index(x,proc_pts[1]-1);
				i = x+(proc_pts[1]-1)*proc_pts[0];
				T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(T[left(i)]+T[right(i)])*pow(deltas[1],2)+ \
				(T[down(i)]+recv_north[x])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
			}
		}
		if(got_north==1 && got_east==1) {
			//i = index(proc_pts[0]-1,proc_pts[1]-1);
			i = proc_pts[0]-1+(proc_pts[1]-1)*proc_pts[0];
			T[i] = (-1*v[i]*pow((deltas[0]*deltas[1]),2)+(recv_east[proc_pts[1]-1]+T[left(i)])*pow(deltas[1],2)+ \
				(recv_north[proc_pts[0]-1]+T[down(i)])*pow(deltas[0],2))/(2*pow(deltas[0],2)+2*pow(deltas[1],2));
		}


		/*
		 * Check the status of your send and receive requests. 
		 */
		if(mycoord[1]%2==0) {
			if(bound_south==0) {
				MPI_Waitall(2, &req[0], &stati[0]);
				got_south = 1;
			}
			if(bound_north==0) {
				MPI_Waitall(2, &req[2], &stati[2]);
				got_north = 1;
			}
		}
		else {
			if(bound_north==0) {
				MPI_Waitall(2, &req[2], &stati[2]);
				got_north = 1;
			}
			if(bound_south==0) {
				MPI_Waitall(2, &req[0], &stati[0]);
				got_south = 1;
			}
			
		}
		if(mycoord[0]%2==0) {
			if(bound_east==0) {
				MPI_Waitall(2, &req[4], &stati[4]);
				got_east = 1;
			}
			if(bound_west==0) {
				MPI_Waitall(2, &req[6], &stati[6]);
				got_west = 1;
			}
		}
		else {
			if(bound_west==0) {
				MPI_Waitall(2, &req[6], &stati[6]);
				got_west = 1;
			}
			if(bound_east==0) {
				MPI_Waitall(2, &req[4], &stati[4]);
				got_east = 1;
			}
		}

		if(will_break==1)
			break;


		if(count%1000==0) {
			memcpy(test_buffer, T, proc_size*sizeof(double));
		}
		if(count%1000==1 && verify_error(T, test_buffer, proc_size)==1) {
			will_break = 1;
		}

		double err=-1;
		if(count%1000==1) {
			err = get_error(T, test_buffer, proc_size);
			printf("(%d): iteration %d... \t%lf\n", myrank, count, err);
		}	

		//printf("SANITY CHECK: \n");
		//printf("remember_south=%d, rem_north=%d, rem_east=%d, rem_west=%d\n", remember_south, remember_north, remember_east, remember_west);


		count++;

		got_east  = 0;
		got_west  = 0;
		got_south = 0;
		got_north = 0;
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		
	}
	printf("(%d): EXITED\n", myrank);
	fflush(stdout);
	

//		if(myrank==0) {
			printf("\n\n(%d):  Printing T(v)\n", myrank);
			for(int i=0; i<proc_size; i++) {
				if(i%proc_pts[0]==0)
					printf("...\n");
				printf("%lf(%lf), ", T[i], v[i]);
			}		
//		}



/*	if(myrank==0) {
		printf("\nPrinting V\n");
		for(int i=0; i<proc_size; i++) {
			if(i%proc_pts[0]==0)
				printf("...\n");
			printf("%lf\n", v[i]);
		}
		
	}
*/	


/*	sleep(2*myrank);
	int count_y=0;
	printf("\n\n\n#######   P %d; (%d, %d)   ########\n", myrank, mycoord[0], mycoord[1]);
	for(int i=0; i<proc_size; i++) {
		if(i%proc_pts[0]==0)
			printf("\n(%d)...", count_y++);
		printf("%lf(%d), ", T[i], i%proc_pts[0]);
	}
	printf("\n\n\n####### END OF  P %d; (%d, %d)   ########\n", myrank, mycoord[0], mycoord[1]);
*/

	double Xmin = (ranges[0]/dims_procs[0])*mycoord[0];
	double Ymin = (ranges[1]/dims_procs[1])*mycoord[1];
	double Xmax = Xmin+(ranges[0]/dims_procs[0]);	
	double Ymax = Ymin+(ranges[1]/dims_procs[1]);
	VTK_out(proc_pts[0], proc_pts[1], &Xmin, &Xmax, &Ymin, &Ymax, T, myrank);

	free(v);
	free(T);
	free(test_buffer);
	free(send_south);
	free(send_north);
	free(send_east);
	free(send_west);
	free(recv_south);
	free(recv_north);
	free(recv_east);
	free(recv_west);


	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
	return 0;
}
