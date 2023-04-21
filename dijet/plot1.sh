#! /usr/bin/env bash

for i in GBW GBWr GBWS GBWSr BGK  BGKr BGKS BGKSr KS
do
	rm ${i} -r
	./run.sh prepare input${i} ${i}
done
for i in GBW GBWr GBWS GBWSr BGK  BGKr BGKS BGKSr KS
do
	cd ${i}
	./optimize.sh
	cd ..
done



