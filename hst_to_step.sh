#! /usr/bin/env bash 

#for i0 in dijet-HERA dijet-HERA-scale-2 dijet-HERA-scale-05
#do
mkdir dijet-HERA-plot
for i1 in GBW BGK GBWS BGKS
do
	for i2 in r kt
	do 
		mkdir dijet-HERA-plot/${i1}${i2}
		for i in {1..8}
		do
			./a.out dijet-HERA/${i1}${i2}/fig16${i}.hst \
				dijet-HERA-scale-2/${i1}${i2}/fig16${i}.hst  \
				dijet-HERA-scale-05/${i1}${i2}/fig16${i}.hst \
				dijet-HERA-plot/${i1}${i2}/fig16${i}.hst 
		done
	done 
done
