#! /usr/bin/env sh

SAT=/home/tomoki/Saturation-Model/saturation-ver3

for i in  ./Run2/fixa-bjorx \
	./Run2/runa-bjorx-4 \
	./Run2/fixa-bjorx-BGK \
	./Run2/GBW ./Run2/GBW-Massive \
	 ./Run2/BGK
do
	echo ${i}
 	export DIR=${i}
	make plot
	${DIR}/grid
	#${DIR}/plot -Q2 100 -x 4
	${DIR}/plot -Q2 100 -x 6
	${DIR}/plot -Q2 100 -x 2
done

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW.png" \
	 -m ./Run2/GBW/gluon-grid.txt ./Run2/fixa-bjorx/gluon-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run.png"  \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/runa-bjorx-4/gluon-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK.png"  \
	-m ./Run2/BGK/gluon-grid.txt ./Run2/fixa-bjorx-BGK/gluon-grid.txt #./Run2/fixa-bjorx-BGK2/gluon-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-contour.png"  \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/fixa-bjorx/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run-contour.png" \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/runa-bjorx-4/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK-contour.png" -c \
	-m ./Run2/BGK/gluon-grid.txt ./Run2/fixa-bjorx-BGK/gluon-grid.txt #./Run2/fixa-bjorx-BGK2/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW.png"  -d \
	-m ./Run2/GBW/dipole-grid.txt ./Run2/fixa-bjorx/dipole-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW-run.png" -d \
	-m ./Run2/GBW/dipole-grid.txt ./Run2/runa-bjorx-4/dipole-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/dipole-BGK.png"  -d \
	-m ./Run2/BGK/dipole-grid.txt  ./Run2/fixa-bjorx-BGK/dipole-grid.txt #./Run2/fixa-bjorx-BGK2/dipole-grid.txt


exit 0

