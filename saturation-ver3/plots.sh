#! /usr/bin/env sh

SAT=/home/tomoki/Saturation-Model/saturation-ver3

#for i in  ./Run3/fixa-bjorx-BGK \
#	./Run3/fixa-bjorx \
#	./Run3/runa-bjorx-4 \
#	./Run3/GBW \
#	 ./Run3/BGK
#do
#	echo ${i}
#	export DIR=${i}
#	make plot
#	
#done

#exit

#for i in  ./Run3/fixa-bjorx-BGK \
#	 ./Run3/BGK \
#	./Run3/fixa-bjorx \
#	./Run3/runa-bjorx-4 \
#	./Run3/GBW 
for i  in ./Run3/GBW \
	./Run3/BGK
do
	export DIR=${i}
	make plot
	${DIR}/F2
	${DIR}/grid
	${DIR}/plot -Q2 100 -x 6
	${DIR}/plot -Q2 100 -x 2
done

#./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW.png" \
#	 -m ./Run3/GBW/gluon-grid.txt ./Run3/fixa-bjorx/gluon-grid.txt 

#./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run.png"  \
#	-m ./Run3/GBW/gluon-grid.txt ./Run3/runa-bjorx-4/gluon-grid.txt 

#./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK.png"  \
#	-m ./Run3/BGK/gluon-grid.txt ./Run3/fixa-bjorx-BGK/gluon-grid.txt #./Run3/fixa-bjorx-BGK2/gluon-grid.txt

#./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-contour.png"  \
#	-m ./Run3/GBW/gluon-grid.txt ./Run3/fixa-bjorx/gluon-grid.txt -c

#./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run-contour.png" \
#	-m ./Run3/GBW/gluon-grid.txt ./Run3/runa-bjorx-4/gluon-grid.txt -c

#./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK-contour.png" -c \
#	-m ./Run3/BGK/gluon-grid.txt ./Run3/fixa-bjorx-BGK/gluon-grid.txt #./Run3/fixa-bjorx-BGK2/gluon-grid.txt -c

#./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW.png"  -d \
#	-m ./Run3/GBW/dipole-grid.txt ./Run3/fixa-bjorx/dipole-grid.txt 

#./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW-run.png" -d \
#	-m ./Run3/GBW/dipole-grid.txt ./Run3/runa-bjorx-4/dipole-grid.txt

#./Plotting/surface.py -s "./kt-formula-report1/dipole-BGK.png"  -d \
#	-m ./Run3/BGK/dipole-grid.txt  ./Run3/fixa-bjorx-BGK/dipole-grid.txt #./Run3/fixa-bjorx-BGK2/dipole-grid.txt


exit 0

