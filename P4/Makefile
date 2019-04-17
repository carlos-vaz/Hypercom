CC = nvcc

CFLAGS = -arch=sm_60

TGT = main

% : %.c
	$(CC) $< -o $@ $(CFLAGS)

all : $(TGT)

clean:
	rm $(TGT)
