# Project 2: 2D Poisson Solver using MPI
Carlos Vazquez Gomez
March 15, 2019



## Description of the files: 

### main.c	
2D Poisson Solver using MPI. Reads .txt files from Input_Data directory, 
and outputs .vtk files into VTK_Data directory. 
Usage:	
```
mpirun -np [P] FILENAME [Py] [Px]
```

### writedata.c
Writes data in the form x*e^y to .txt file in Input_Data directory. 
Usage:	
```
./writedata [Y Range] [X Range] [# points y] [# points x]
```

### vtk.c 
Defines extern function linked into 'main' to create .vtk files. 
Usage:	N/A (linked into 'main')
	

### Makefile
Compiles and links 'main.c', 'writedata.c', and 'vtk.c' into two executables
'main' and 'writedata'. 
Usage:	
```
make clean
make all
```

### Slurm/
	Directory with slurm scripts for weak scaling, convergence, grid convergence, 
	and heterogenous grid tests. The output of each script is a .out file which
	can be parsed by awk to retrieve timing and cycle counts. 
	Usage (example):
```	
sbatch Slurm/slurm_weak_analysis
grep weak_analysis.out -e "time" | awk -F"= " '{ print $2 }'
grep weak_analysis.out -e "cycles" | awk -F"= " '{ print $2 }'
```
