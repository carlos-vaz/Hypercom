CC = gcc

TGT = main

CFLAGS = -lm -std=c99

LFLAGS = -fopenmp

% : %.c
	$(CC) $< vtk.c -o $@ $(CFLAGS) $(LFLAGS)

all: $(TGT)

clean:
	rm $(TGT)