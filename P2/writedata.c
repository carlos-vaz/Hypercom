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
	printf("x_dim=%d\ny_dim=%d\nWriting data to data_%d_%d.txt\n",x_dim, y_dim, x_dim, y_dim);
	int fd = open("./data_%d.txt", O_WRONLY | O_CREAT, 0644, 1);
	char * test = "hello\n";
	write(fd, test, 6);
	close(fd);
}


