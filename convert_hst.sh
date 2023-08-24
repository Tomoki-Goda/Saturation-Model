#! /usr/bin/env bash


for i in GBWS GBW BGKS BGKS
do 
	for j in r kt
	do
		for k in {1..8}
		do
			./a.out ./dijet-HERA/${i}${j}/fig16${k}.hst \
				./dijet-HERA-scale-2/${i}${j}/fig16${k}.hst \
				./dijet-HERA-scale-05/${i}${j}/fig16${k}.hst \
				./dijet-HERA-plot/${i}${j}/fig16${k}.hst
		done
	done
done



