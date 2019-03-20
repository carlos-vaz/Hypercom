### Numerical Integration with Posix threads

Example: The following commands compile main and invoke it to integrate 4/(1+x^2) from 0 to 1 using 4000000 points and 8 threads. 
```
$ make
[BUILDING main...]
$ ./main 4000000 8
Running 8 threads on 4000000 points
FINAL VALUE: 3.141593
ELAPSED TIME: 0.003656
$ 
```

You can also run a strong scaling analysis script, which tests on 4
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
