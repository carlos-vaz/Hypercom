#include "stdlib.h"
#include "stdio.h"
#include "string.h"

int xdim, ydim;

int main(int argc char **argv) {
	if(argc!=3) {
		printf("Usage:\n\t./main [xdim] [ydim]\n");
	}
	xdim = atoi(argv[1]);
	ydim = atoi(argv[2]);
}
