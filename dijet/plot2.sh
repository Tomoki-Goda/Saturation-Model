#! /usr/bin/env bash

#it is recommended one switch from read_event_file to create_event_file

for i in GBW GBWr GBWS GBWSr BGK  BGKr BGKS BGKSr KS
do
	cd ${i}
	parallel ./main.out :::  seed=2023 seed=5432 seed=1357 seed=2468 
	#./main.out   seed=1234 
	rm eventfile.dat Phi*
	../run.sh merge raw*
	../read_event_file.out eventfile.dat
	echo ${i} ' done'
	cd ..
done



