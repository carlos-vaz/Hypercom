## Numerical Integration with Posix threads

### Build
```
$ make
[BUILDING main...]
$
```

### Run
Example:  
Integrate 4/(1+x^2) from 0 to 1 using (defined internally in main.c) 4000000 points and 8 threads:
$ ./main 4000000 8
Running 8 threads on 4000000 points
FINAL VALUE: 3.141593
ELAPSED TIME: 0.003656
$


### Strong Scaling Analysis
You can also run a strong scaling analysis script:
```
$ bash scaling.sh 
Testing on 400000000 points with 1, 2, 4, 8, 16, 32, 64 threads...
Progress:   |################################################# |

Timings:
____________________
7.946787
4.803657
2.523130
5.185251
1.196787
1.107933
0.918316

Values:
____________________
3.1415926561
3.1415926561
3.1415926561
3.1415926561
3.1415926561
3.1415926561
3.1415926561

$
```
