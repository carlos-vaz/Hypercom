CC = mpicc
CFLAGS = -std=c99
LFLAGS = -lm
OBJ =	main \
	native \
	slow

comp = 	 @module load openmpi; \
	echo '[Compiling $@...]'; \
	$1;

% : %.c
	$(call comp, $(CC) $(CFLAGS) $(LFLAGS) $< -o $@)

all: $(OBJ)

clean: 
	rm $(OBJ)
