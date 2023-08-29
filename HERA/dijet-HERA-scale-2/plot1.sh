#! /usr/bin/env bash

for i in BGKkt  BGKSkt
do
	rm ${i} -r -f
	./run.sh prepare input${i} ${i}
done
for i in BGKkt  BGKSkt
do
	cd ${i}
	./optimize.sh Nparallel=6
	cd ..
done



