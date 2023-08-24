#! /usr/bin/env bash

for i in GBWkt GBWr GBWSkt GBWSr BGKkt  BGKr BGKSkt BGKSr KS rcBK
do
	cd ${i}
  	#parallel ./main.out :::  seed=2023 seed=5432 seed=1357 seed=2468
	parallel ./main.out :::  seed=1904 seed=794 seed=1192 seed=1702
	#./main.out   seed=1234 
	rm -f eventfile.dat *.hst
	#../run.sh merge raw*
	../create_eventfile.out lhef,pb raw*
	echo ${i} ' done'
	cd ..
done



