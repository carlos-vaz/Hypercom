#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

int nt, np;

int main(int argc, char *argv[]) {
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
	return 0;
}
