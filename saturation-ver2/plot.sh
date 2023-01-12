#!/usr/bin/env bash

#for dir in ../Saturation-Model2207/Run2207*/*/M* # ../Results/BGK/M* ../Results/GBWS/M* ../Results/BGKS/M*
for dir in ../Run/*/M*
do
	for x in 2 4 6
	do
		q2=0
#		for q2 in 5 50 100 500
#		do		
#			${dir}/dipole -in ${dir}/result.txt -out ${dir}/dipole-${q2}-${x}.txt -x ${x} -Q2 ${q2}
			${dir}/tmd-gluon -in ${dir}/result.txt -out ${dir}/gluon-${q2}-${x}.txt -x ${x} -Q2 ${q2}
			${dir}/ww-gluon -in ${dir}/result.txt -out ${dir}/ww-gluon-${q2}-${x}.txt -x ${x} -Q2 ${q2}
#			#${dir}/Integrand -in ${dir}/result.txt -out ${dir}/integrand-${q2}-${x}.txt -x ${x} -Q2 ${q2}
		#done
		
#		${dir}/F2 -in ${dir}/result.txt -out ${dir}/F2-${x}.txt -x ${x}
#		${dir}/F2-slope -in ${dir}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x}
			
	done
	q2=0
#	for q2 in 5 50 500
#	do	
		for k in 1
		do
			${dir}/tmd-gluon-x -in ${dir}/result.txt -out ${dir}/gluon-x-${q2}-${k}.txt -Q2 ${q2} -k ${k}
			${dir}/ww-gluon-x -in ${dir}/result.txt -out ${dir}/ww-gluon-x-${q2}-${k}.txt -Q2 ${q2} -k ${k}
		done
		#${dir}/critical -in ${dir}/result.txt -out ${dir}/critical-${q2}.txt -Q2 ${q2}
		${dir}/tmd-critical -in ${dir}/result.txt -out ${dir}/tmd-critical-${q2}.txt -Q2 ${q2}
		#${dir}/ww-critical -in ${dir}/result.txt -out ${dir}/ww-critical-${q2}.txt -Q2 ${q2}
		#echo "hello"
#	done
#	${dir}/fcn -in ${dir}/result.txt -out ${dir}/fcn.txt
#	${dir}/export-data -in ${dir}/result.txt -out ${dir}/data.txt	
done

	
