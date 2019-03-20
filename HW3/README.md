## Numerical Integration with Posix threads

### Build
```
$ make
[BUILDING main...]
$
```

### Run
Example:  
Integrate 4/(1+x^2) from 0 to 1 (defined internally in main.c) using 4000000 points and 8 threads:
```
$ ./main 4000000 8
Running 8 threads on 4000000 points
FINAL VALUE: 3.141593
ELAPSED TIME: 0.003656
$
```


### Strong Scaling Analysis
You can also run a strong scaling analysis script:
```
$ bash scaling.sh 
Testing on 400000000 points with 1, 2, 4, 8, 16, 32, 64 threads...
Progress:   |################################################# |

Timings:
____________________
ELAPSED TIME (1 threads): 8.528089
ELAPSED TIME (2 threads): 4.873009
ELAPSED TIME (4 threads): 2.846566
ELAPSED TIME (8 threads): 1.717208
ELAPSED TIME (16 threads): 1.197903
ELAPSED TIME (32 threads): 1.026930
ELAPSED TIME (64 threads): 0.967940

Values:
____________________
FINAL VALUE (1 threads): 3.1415926561
FINAL VALUE (2 threads): 3.1415926561
FINAL VALUE (4 threads): 3.1415926561
FINAL VALUE (8 threads): 3.1415926561
FINAL VALUE (16 threads): 3.1415926561
FINAL VALUE (32 threads): 3.1415926561
FINAL VALUE (64 threads): 3.1415926561

$
```
