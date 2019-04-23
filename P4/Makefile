CC = nvcc

CFLAGS = -arch=sm_60 \
	 -rdc=true

TGT = main \
  	coop

% : %.cu
	$(CC) $< vtk.c -o $@ $(CFLAGS)

all : $(TGT)

clean:
	rm $(TGT)
