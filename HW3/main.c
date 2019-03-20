#define _POSIX_C_SOURCE 199309L

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

#define START	0
#define END	1

double * buf, sum;
int nt;
long np;
pthread_mutex_t sum_mutex = PTHREAD_MUTEX_INITIALIZER;


/* 
 * Routine that each thread executes to compute its share
 * of the integral
 */ 
void *thread_routine(void *ID) {
	//printf("Thread %d checking in\n", *(int*)ID);
	double mysum = 0;
	double del_x = (double)(END-START)/np;
	for(int i=(*(int*)ID)*np/nt; i<(*(int*)ID+1)*np/nt; i++) {
		mysum += buf[i]*del_x;
	}
	pthread_mutex_lock(&sum_mutex);
	sum += mysum;
	pthread_mutex_unlock(&sum_mutex);	
}


/*
 * Interchangeable function that defines curve to integrate
 */
double function(double x) {
	return 4/(1+pow(x,2));
}

int main(int argc, char *argv[]) {
	/* 
	 * Parse arguments
	 */
	if(argc!=3) {
		printf("Usage:\n\t./main [# points] [# threads]\n");
		exit(0);
	}
	np = strtol(argv[1], NULL, 10);
	nt = atoi(argv[2]);
	if(np%nt!=0) {
		printf("Error: # threads must divide # points\n");
		exit(1);
	}
	printf("Running %d threads on %ln points\n", nt, np);

	/* 
	 * Build function
	 */
	buf = (double*)malloc(np*sizeof(double));
	double x;
	for(int i=0; i<np; i++) {
		buf[i] = function( (double)i/np );
	}
	
	/* 
	 * Dispatch threads
	 * Start timer
	 */
	struct timespec start, end; 
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC, &start);
	pthread_t *threads = (pthread_t*)malloc(nt*sizeof(pthread_t));
	int *ID = malloc(nt*sizeof(int));
	
	for(int i=0; i<nt; i++) {
		ID[i] = i;
		pthread_create(&threads[i], NULL, thread_routine, &ID[i]);
	}
	for(int i=0; i<nt; i++) {
		pthread_join(threads[i], NULL);
	}
	printf("FINAL VALUE: %.10f\n", sum);

	/* 
	 * Display elapsed time
	 */
	clock_gettime(CLOCK_MONOTONIC, &end);
	elapsed = end.tv_sec - start.tv_sec;
	elapsed += (end.tv_nsec - start.tv_nsec)/1000000000.0;
	printf("ELAPSED TIME: %lf\n", elapsed);

	return 0;
}
