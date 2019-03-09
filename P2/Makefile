CC = mpicc
CFLAGS = -std=c99
LFLAGS = -lm
OBJ =	main \
	writedata

comp =  \
	@echo '[Compiling $@...]'; \
	$1;
link = \
	@echo '[Linking $@...]'; \
	$1;

writedata : writedata.c
	$(call comp, gcc writedata.c -o writedata $(LFLAGS))

main : main.c vtk.c
	$(call link, $(CC) $(CFLAGS) $< vtk.c -o $@ $(LFLAGS))

all: $(OBJ)

clean: 
	rm $(OBJ)
