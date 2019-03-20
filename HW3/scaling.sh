#!/usr/bin/env bash

echo > log.out
for i in {1 2 4 8 16 32 64}
do
	./main 40000000 $i >> log.out
	let "n = $( (echo 'l(8)/l(2)' | bc -l) | cut -c1-1)"
	echo $n
done

echo "Timings:
____________________
"
grep log.out -e "ELAPSED" | awk -F": " '{print $2}' 

echo "
Values:
____________________
"
grep log.out -e "FINAL" | awk -F": " '{print $2}' 

