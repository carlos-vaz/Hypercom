#!/usr/bin/env bash

echo > log.out
echo "Testing on 400000000 points with 1, 2, 4, 8, 16, 32, 64 threads..."
printf "\xd"
printf "Progress:   |"
printf " %.0s" $(seq 1 $((7*7)))
printf "|"
for i in 1 2 4 8 16 32 64
do
	./main 400000000 $i >> log.out
	let "n = $( (echo 'l('$i')/l(2)' | bc -l) | cut -c1-1) + 1"
	printf "\xd"
	printf "Progress:   |"
	printf "#%.0s" $(seq 1 $((n*7)) )
	printf " %.0s" $(seq 1 $(( (7-n)*7 )) )
	printf "|"
done

echo "

Timings:
____________________"
grep log.out -e "ELAPSED" #| awk -F": " '{print $2}' 

echo "
Values:
____________________"
grep log.out -e "FINAL" #| awk -F": " '{print $2}' 
echo 

