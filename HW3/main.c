#include <stdio.h>
#include <pthread.h>

int nt, np;

int main(int argc, char *argv[]) {
	if(argc!=3) {
		printf("Usage:\n\t./main [# points] [# threads]\n");
		exit(0);
	}
	np = atoi(argv[1]);
	nt = atoi(argv[2]);
	printf("Running %d threads on %d points\n", nt, np);
	return 0;
}
