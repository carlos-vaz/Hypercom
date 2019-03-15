Project 2: 2D Poisson Solver using MPI
Carlos Vazquez Gomez
March 15, 2019



Description of the files: 

main.c	
	2D Poisson Solver using MPI. Outputs .vtk files into VTK_Data directory. 
	Usage:	mpirun -np X FILENAME Py Px

writedata.c
	Writes data in the form x*e^y to .txt file in Input_Data directory. 
	Usage:	./writedata [Y Range] [X Range] [# points y] [# points x]

vtk.c 
	Defines extern function linked into 'main' to create .vtk files. 
	Usage:	N/A (linked into 'main')
	

Makefile
	Compiles and links 'main.c', 'writedata.c', and 'vtk.c' into two executables
	'main' and 'writedata'. 
	Usage:	make clean; make all;

crc_slurm
	Slurm script that allocates 100 processes on H2P cluster and loops through 
	the weak scaling analysis tests. The output is Toro.out, which can be parsed
	by awk to retrieve timing and cycle counts. 
	Usage:	sbatch crc_slurm
		grep Toro.out -e "time" | awk -F"= " '{ print $2 }'
		grep Toro.out -e "cycles" | awk -F"= " '{ print $2 }'

