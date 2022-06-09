#! /usr/bin/env bash


gcc ./Utilities/Auto-Plot-DP.c -o ./Auto-DP -lm
gcc ./Utilities/Auto-Plot-F2.c -o ./Auto-F2 -lm

dir="NewRun"
for x in 3 5
do 
	for q2 in 50 500 #10 50 100 500
	do
		./Auto-DP ${q2} ${x} ./${dir}/*/Mass*
	done
	./Auto-F2 ${x} ./${dir}/*/Mass*
done 

rm ./Auto-F2
rm ./Auto-DP
