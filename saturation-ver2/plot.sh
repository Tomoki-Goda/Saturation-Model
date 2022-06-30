#!/usr/bin/env bash

for dir in ./GBW/M* ./BGKStar/M* ./GBWSStar/M* ./BGKSStar/M*
do
	for x in 3 4 5
	do
		for q2 in 10 50 500
		do		
			${dir}/dipole -in ${dir}/result.txt -out ${dir}/dipole-${q2}-${x}.txt -x ${x} -Q2 ${q2}
		done
		
		${dir}/F2 -in ${dir}/result.txt -out ${dir}/F2-${x}.txt -x ${x}
		${dir}/F2-slope -in ${dir}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x}
			
	done
	for q2 in 10 50 500
		do		
			${dir}/critical -in ${dir}/result.txt -out ${dir}/critical-${q2}.txt -Q2 ${q2}
		done
	
done
