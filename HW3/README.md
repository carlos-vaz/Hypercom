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
Progress:   |#####################                            |
```

and eventually: 
```
$ bash scaling.sh 
Running:  |################################### |

Timings:
____________________
0.475256
0.259443
0.132217
0.066913
0.061596
0.047071
0.050738

Values:
____________________
3.1415926786
3.1415926786
3.1415926786
3.1415926786
3.1415926786
3.1415926786
3.1415926786

$
```
