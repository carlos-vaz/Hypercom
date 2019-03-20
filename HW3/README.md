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

You can also run a strong scaling analysis script. 
```

```
