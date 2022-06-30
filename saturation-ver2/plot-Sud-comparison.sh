#!/usr/bin/env bash

export DIR

for i in ./BGK-BGKS-Star/M* ./GBW-GBWS-Star/M*;
do
       	DIR=${i};
	make;
done


for q2 in 50 500;
do
	for m in 0.0 0.0196;
	do
		dir=./BGK-BGKS-Star/Mass${m}-Qup${q2}-Model1-Sud0-rfix0;
		dirS=./BGKSStar/Mass${m}-Qup${q2}-Model3-Sud1-rfix0;
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



for q2 in 50 500;
do
	for m in 0.0 0.0196;
	do
		dir=./GBW-GBWS-Star/Mass${m}-Qup${q2}-Model0-Sud0-rfix0;
		dirS=./GBWSStar/Mass${m}-Qup${q2}-Model22-Sud1-rfix0;
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
