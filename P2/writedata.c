#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>

/* 
 * Writes exponential into grid
 */

int main(int argc, char* argv[]) {
	if(argc!=3) {
		printf("Not correct # of args\n");
		exit(0);
	}

	int x_dim = atoi(argv[1]);
	int y_dim = atoi(argv[2]);
	char * filename = malloc(20);
	filename = sprintf(filename, "data_%dx%d.txt", x_dim, y_dim);
	printf("x_dim=%d\ny_dim=%d\nWriting data to %s",x_dim, y_dim, filename);
	int fd = fopen(filename, O_WRONLY | O_CREAT, 0644, 1);
	char * test = "hello\n";
	write(fd, test, 6);
	close(fd);
}


