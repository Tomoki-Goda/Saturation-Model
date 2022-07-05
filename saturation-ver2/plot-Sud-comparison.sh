#!/usr/bin/env bash

export DIR

for i in ../Results/BGK-BGKS/M* ../Results/GBW-GBWS/M*;
do
       	DIR=${i};
	make;
done


for q2 in 10 50 650;
do
	for m in 0.0;
	do
		dir=../Results/BGK-BGKS/Mass${m}-Qup${q2}-Model1-Sud0-rfix0;
		dirS=../Results/BGKS/Mass${m}-Qup${q2}-Model3-Sud1-rfix0;
		for x in 3 5;
		do
			for Q2 in 10 50 500;
			do
				${dir}/dipole -in ${dirS}/result.txt -out ${dir}/dipole-${Q2}-${x}.txt -Q2 ${Q2} -x ${x};
			done;
	
			${dir}/F2 -in ${dirS}/result.txt -out ${dir}/F2-${x}.txt -x ${x};
			${dir}/F2-slope -in ${dirS}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x};
		done;

		for Q2 in 10 50 500 ;
		do
			${dir}/critical -in ${dirS}/result.txt -out ${dir}/critical-${Q2}.txt -Q2 ${Q2};
		done;
	done;
done



for q2 in 10 50 650;
do
	for m in 0.0;
	do
		dir=../Results/GBW-GBWS/Mass${m}-Qup${q2}-Model0-Sud0-rfix0;
		dirS=../Results/GBWS/Mass${m}-Qup${q2}-Model22-Sud1-rfix0;
		for x in 3 5;
		do
			for Q2 in 10 50 500;
			do
				${dir}/dipole -in ${dirS}/result.txt -out ${dir}/dipole-${Q2}-${x}.txt -Q2 ${Q2} -x ${x};
			done;
	
			${dir}/F2 -in ${dirS}/result.txt -out ${dir}/F2-${x}.txt -x ${x};
			${dir}/F2-slope -in ${dirS}/result.txt -out ${dir}/F2-slope-${x}.txt -x ${x};
		done;

		for Q2 in 10 50 500 ;
		do
			${dir}/critical -in ${dirS}/result.txt -out ${dir}/critical-${Q2}.txt -Q2 ${Q2};
		done;
	done;
done
