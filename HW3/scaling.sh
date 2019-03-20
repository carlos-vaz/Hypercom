#!/usr/bin/env bash

echo > log.out
for i in {1 2 4 8 16 32 64}
do
	./main 40000000 $i >> log.out
	let "n = $( (echo 'l('$i')/l(2)' | bc -l) | cut -c1-1)"
	printf "\xd"
	printf "|"
	printf "#%.0s" $(seq 1 $((n*3)) )
	printf " %.0s" $(seq 1 $( (6-n)*3 ) )
	printf "|"
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

