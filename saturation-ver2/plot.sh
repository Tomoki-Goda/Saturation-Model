#!/usr/bin/env bash

#for dir in ../Saturation-Model2207/Run2207*/*/M* # ../Results/BGK/M* ../Results/GBWS/M* ../Results/BGKS/M*
for dir in ../Run2607-Fejer/*/M* 
do
	for x in 3 5 7
	do
		for q2 in 100 650
		do		
			${dir}/dipole -in ${dir}/result.txt -out ${dir}/dipole-${q2}-${x}.txt -x ${x} -Q2 ${q2}
		done
		
		${dir}/F2 -in ${dir}/result.txt -out ${dir}/F2-${x}.txt -x ${x}
		${dir}/F2-slope -in ${dir}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x}
			
	done
	for q2 in 100 650
		do		
			${dir}/critical -in ${dir}/result.txt -out ${dir}/critical-${q2}.txt -Q2 ${q2}
		done
	
done
