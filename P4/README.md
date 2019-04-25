### Solving 2D Poisson with CUDA

To run on h2p cluster:
```
./prepare
module load cuda/9.0
make
./main 500 500
'''

## Files
# main.cu
Solves grid of any size by launching multiple kernels in a loop
and synchronizing with cudaDeviceSynchronize.

# coop.cu
Solves up to 56,000 point grids (e.g. 200x280) by synchronizing 
the whole grid using a Cooperative Kernel. Much faster than main.cu. 

