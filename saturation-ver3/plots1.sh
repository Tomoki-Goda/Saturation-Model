#! /usr/bin/env sh

SAT=/home/tomoki/saturation-ver2/Run2

for i in  ${SAT}/*/Mass* 
do
	echo ${i}
	export DIR=${i}
	make plot
done


for i in  ${SAT}/*/Mass* 
do
	export DIR=${i}
	#make plot
	${DIR}/F2
	${DIR}/grid
	${DIR}/plot -Q2 100 -x 6
	${DIR}/plot -Q2 100 -x 2
done

exit 0

