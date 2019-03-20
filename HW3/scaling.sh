#!/usr/bin/env bash

for i in 1 2 4 8 16 32 64
do
	./main 4000000 $i >> log.out
done

echo "Timings:\n____________________"
grep log.out -e "ELAPSED" | awk -F": " '{print $2}' 

echo "\nValues:\n____________________"
grep log.out -e "FINAL" | awk -F": " '{print $2}' 

