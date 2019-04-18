CC = nvcc

CFLAGS = -arch=sm_60

TGT = main

% : %.cu
	$(CC) $< -o $@ $(CFLAGS)

all : $(TGT)

clean:
	rm $(TGT)
