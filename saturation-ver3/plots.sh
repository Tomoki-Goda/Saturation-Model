#! /usr/bin/env sh

SAT=/home/tomoki/Saturation-Model/saturation-ver3

for i in  ./Run3/fixa-bjorx-BGK \
	./Run3/fixa-bjorx \
	./Run3/runa-bjorx-4 \
	./Run3/GBW \
	 ./Run3/BGK
do
	echo ${i}
	export DIR=${i}
	make plot
done


for i in  ./Run3/fixa-bjorx-BGK \
	 ./Run3/BGK \
	./Run3/fixa-bjorx \
	./Run3/runa-bjorx-4 \
	./Run3/GBW  
do
	export DIR=${i}
	#make plot
	${DIR}/F2
	${DIR}/grid
	${DIR}/plot -Q2 100 -x 6
	${DIR}/plot -Q2 100 -x 2
done

exit 0

