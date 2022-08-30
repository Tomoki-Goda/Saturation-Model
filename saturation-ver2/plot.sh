#!/usr/bin/env bash

#for dir in ../Saturation-Model2207/Run2207*/*/M* # ../Results/BGK/M* ../Results/GBWS/M* ../Results/BGKS/M*
for dir in ../Run2808/*/M*

do
	for x in 2 4 6
	do
		for q2 in 5 100 650
		do		
			${dir}/dipole -in ${dir}/result.txt -out ${dir}/dipole-${q2}-${x}.txt -x ${x} -Q2 ${q2}
			${dir}/tmd-gluon -in ${dir}/result.txt -out ${dir}/gluon-${q2}-${x}.txt -x ${x} -Q2 ${q2}
			${dir}/Integrand -in ${dir}/result.txt -out ${dir}/integrand-${q2}-${x}.txt -x ${x} -Q2 ${q2}
		done
		
		${dir}/F2 -in ${dir}/result.txt -out ${dir}/F2-${x}.txt -x ${x}
		${dir}/F2-slope -in ${dir}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x}
			
	done
	for q2 in 5 100 650
	do	
		for k in 5 100 650
		do
			${dir}/tmd-gluon-x -in ${dir}/result.txt -out ${dir}/gluon-x-${k}-${q2}.txt -Q2 ${q2} -k ${k}
		done
		${dir}/critical -in ${dir}/result.txt -out ${dir}/critical-${q2}.txt -Q2 ${q2}
		${dir}/tmd-critical -in ${dir}/result.txt -out ${dir}/tmd-critical-${q2}.txt -Q2 ${q2}
		#echo "hello"
	done
	
	${dir}/fcn -in ${dir}/result.txt -out ${dir}/fcn.txt
	${dir}/export-data -in ${dir}/result.txt -out ${dir}/data.txt
		
done

	
