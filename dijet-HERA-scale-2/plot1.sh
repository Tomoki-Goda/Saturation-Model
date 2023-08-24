#! /usr/bin/env bash

for i in GBWkt GBWr GBWSkt GBWSr BGKkt  BGKr BGKSkt BGKSr KS rcBK
do
	rm ${i} -r -f
	./run.sh prepare input${i} ${i}
done
for i in GBWkt GBWr GBWSkt GBWSr BGKkt  BGKr BGKSkt BGKSr KS rcBK
do
	cd ${i}
	./optimize.sh Nparallel=2
	cd ..
done



