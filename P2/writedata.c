#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>

/* 
 * Writes exponential into grid
 */

int main(int argc, char* argv[]) {
	if(argc!=3) {
		printf("ERROR! Format: ./writedata [# points y] [# points x]\n");
		exit(0);
	}

	int x_dim = atoi(argv[2]);
	int y_dim = atoi(argv[1]);
	char * filename = malloc(20);
	sprintf(filename, "data_%dx%d.txt", y_dim, x_dim);
	printf("x_dim=%d\ny_dim=%d\nWriting data to %s\n",x_dim, y_dim, filename);
	int fd = open(filename, O_WRONLY | O_CREAT, 0644, 1);


	// Write dimensions metadata into file
	write(fd, &x_dim, sizeof(int));
	write(fd, &y_dim, sizeof(int));


	// Write data into file
	double e = 2.718281828, val; 
	for(double y=0; y<y_dim; y++)
		for(double x=0; x<x_dim; x++) {
			//val = x*pow(e, y);
			val = y*x_dim+x;
			write(fd, &val, sizeof(double));
		}
	close(fd);
}


