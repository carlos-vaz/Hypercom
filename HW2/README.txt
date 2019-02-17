'main.c' and 'native.c' are my collectives implementation and MPI's collectives
implementations, resp.  'slow.c' is an earlier attempt in which each process
propagates data only to two child processes. 'main.c' is faster as each process
keeps propagating data down while it cans (like MPI scatter).

Each of these programs takes 3 arguments: the number of points total, the 
starting x value, and the ending x value of integration. 
i.e. mpirun -np [# procs] [program] [# points total] [x_min] [x_max]
To change the function (currently 4/(1+x^2)), change the first function
in each .c file. 

Implementation details can be found in 'Report.pdf'. 

To run, do:

module load openmpi
make
mpirun -np 12 main 12000000 0 1
mpirun -np 12 native 12000000 0 1
mpirun -np 12 slow 12000000 0 1