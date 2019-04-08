CC = gcc

TGT = main

CFLAGS = -lm

LFLAGS = -fopenmp

% : %.c
	$(CC) $< vtk.c -o $@ $(CFLAGS) $(LFLAGS)

all: $(TGT)

clean:
	rm $(TGT)