#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>

int nt, np;

/* 
 * Routine that each thread executes to compute its share
 * of the integral
 */ 
void *thread_routine(void *ID) {
	printf("Thread %d checking in\n", *(int*)ID);
	
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
	np = atoi(argv[1]);
	nt = atoi(argv[2]);
	if(np%nt!=0) {
		printf("Error: # threads must divide # points\n");
		exit(1);
	}
	printf("Running %d threads on %d points\n", nt, np);

	/* 
	 * Build function
	 */
	double *buf = malloc(np*sizeof(double));
	double x;
	for(int i=0; i<np; i++) {
		buf[i] = function( (double)i/np );
	}
	
	/* 
	 * Dispatch threads
	 * Start timer
	 */
	time_t start = clock();

	pthread_t *threads = malloc(nt*sizeof(pthread_t));
	int *ID = malloc(nt*sizeof(int));
	
	for(int i=0; i<nt; i++) {
		ID[i] = i;
		pthread_create(&threads[i], NULL, thread_routine, &ID[i]);
	}

	for(int i=0; i<nt; i++) {
		pthread_join(threads[i], NULL);
	}

	return 0;
}
