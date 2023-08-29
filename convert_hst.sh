#! /usr/bin/env bash

mkdir dijet-HERA-plot
for i in BGKS BGK GBWS GBW
do 
	for j in  kt r
	do
		mkdir dijet-HERA-plot/${i}${j}
		for k in {1..8}
		do
			./a.out ./HERA/dijet-HERA/${i}${j}/fig16${k}.hst \
				./HERA/dijet-HERA-scale-2/${i}${j}/fig16${k}.hst \
				./HERA/dijet-HERA-scale-05/${i}${j}/fig16${k}.hst \
				./dijet-HERA-plot/${i}${j}/fig16${k}.hst
		done
	done
done



