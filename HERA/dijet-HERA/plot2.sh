#! /usr/bin/env bash

for i in GBWkt GBWr BGKkt  BGKr BGKSkt
do
	cd ${i}
 	parallel -j 6 ./main.out :::  seed=2023 seed=5432 seed=1357 seed=2468 seed=1904 seed=794 seed=1192 seed=1702 seed=538 seed=710 seed=1994 seed=1936

	#./main.out   seed=1234 
	rm -f eventfile.dat *.hst
	#../run.sh merge raw*
	../create_eventfile.out lhef,pb raw*
	echo ${i} ' done'
	cd ..
done



