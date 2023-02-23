#! /usr/bin/env sh

SAT=/home/tomoki/Saturation-Model/saturation-ver3

for i in  ./Run2/fixa-bjorx ./Run2/fixa-bjorx2 \
	./Run2/runa-bjorx-4 ./Run2/runa-bjorx2-4 \
	./Run2/fixa-bjorx-BGK ./Run2/fixa-bjorx-BGK2 \
	./Run2/GBW ./Run2/GBW-Massive ./Run2/BGK ./Run2/BGK
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
	 -m ./Run2/GBW/gluon-grid.txt ./Run2/fixa-bjorx/gluon-grid.txt ./Run2/fixa-bjorx2/gluon-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run.png"  \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/runa-bjorx-4/gluon-grid.txt ./Run2/runa-bjorx2-4/gluon-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK.png"  \
	-m ./Run2/BGK/gluon-grid.txt ./Run2/fixa-bjorx-BGK/gluon-grid.txt #./Run2/fixa-bjorx-BGK2/gluon-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-contour.png"  \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/fixa-bjorx/gluon-grid.txt ./Run2/fixa-bjorx2/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run-contour.png" \
	-m ./Run2/GBW/gluon-grid.txt ./Run2/runa-bjorx-4/gluon-grid.txt ./Run2/runa-bjorx2-4/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK-contour.png" -c \
	-m ./Run2/BGK/gluon-grid.txt ./Run2/fixa-bjorx-BGK/gluon-grid.txt #./Run2/fixa-bjorx-BGK2/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW.png"  -d \
	-m ./Run2/GBW/dipole-grid.txt ./Run2/fixa-bjorx/dipole-grid.txt ./Run2/fixa-bjorx2/dipole-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW-run.png" -d \
	-m ./Run2/GBW/dipole-grid.txt ./Run2/runa-bjorx-4/dipole-grid.txt ./Run2/runa-bjorx2-4/dipole-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/dipole-BGK.png"  -d \
	-m ./Run2/BGK/dipole-grid.txt  ./Run2/fixa-bjorx-BGK/dipole-grid.txt #./Run2/fixa-bjorx-BGK2/dipole-grid.txt



./plot.py -s "./kt-formula-report1/GBW-gluon.png" -y "linear" \
	-c "b r r b r r b r r" -p "-. -. - -. -. - -. -. -" -l "GBW-r:GBW-kt:GBW-kt-run::::::" \
	${SAT}/Run2/GBW/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx/gluon-2-100.txt \
	${SAT}/Run2/runa-bjorx-4/gluon-2-100.txt \
	${SAT}/Run2/GBW/gluon-6-100.txt \
	${SAT}/Run2/fixa-bjorx/gluon-6-100.txt \
	${SAT}/Run2/runa-bjorx-4/gluon-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-dipole.png" \
	-c "b r r b r r b r r" -p "-. -. - -. -. - -. -. -" -l "GBW-r:GBW-kt:GBW-kt-run::::::" \
	${SAT}/Run2/GBW/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx/dipole-2-100.txt \
	${SAT}/Run2/runa-bjorx-4/dipole-2-100.txt \
	${SAT}/Run2/GBW/dipole-6-100.txt \
	${SAT}/Run2/fixa-bjorx/dipole-6-100.txt \
	${SAT}/Run2/runa-bjorx-4/dipole-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-thresh-gluon.png" -y "linear" \
	-c "b r r b r r b r r" -p "-. -. - -. -. - -. -. -" -l "GBW-r:GBW-kt:GBW-kt-run::::::" \
	${SAT}/Run2/GBW/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx2/gluon-2-100.txt \
	${SAT}/Run2/runa-bjorx2-4/gluon-2-100.txt \
	${SAT}/Run2/GBW/gluon-6-100.txt \
	${SAT}/Run2/fixa-bjorx2/gluon-6-100.txt \
	${SAT}/Run2/runa-bjorx2-4/gluon-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-thresh-dipole.png" \
	-c "b r r b r r b r r" -p "-. -. - -. -. - -. -. -" -l "GBW-r:GBW-kt:GBW-kt-run::::::" \
	${SAT}/Run2/GBW/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx2/dipole-2-100.txt \
	${SAT}/Run2/runa-bjorx2-4/dipole-2-100.txt \
	${SAT}/Run2/GBW/dipole-6-100.txt \
	${SAT}/Run2/fixa-bjorx2/dipole-6-100.txt \
	${SAT}/Run2/runa-bjorx2-4/dipole-6-100.txt 
	
./plot.py -s "./kt-formula-report1/BGK-gluon.png" -y "linear" \
	-c "b r b r b r" -p "-. - -. - -. -" -l "BGK-r:BGK-kt::::::" \
	${SAT}/Run2/BGK/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/gluon-2-100.txt \
	${SAT}/Run2/BGK/gluon-6-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/gluon-6-100.txt 
./plot.py -s "./kt-formula-report1/BGK-dipole.png" \
	 -c "b r b r b r" -p "-. - -. - -. -" -l "BGK-r:BGK-kt::::::" \
	${SAT}/Run2/BGK/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/dipole-2-100.txt \
	${SAT}/Run2/BGK/dipole-6-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/dipole-6-100.txt 
	
	
	
#for x in 2 4 6
#do
#./plot.py -s "./kt-formula-report1/GBW-gluon-${x}-100.png" -y "linear" -c "b r r" -p "-. -. -" -l "GBW-r:GBW-kt:GBW-kt-run" \
#	${SAT}/Run2/GBW/gluon-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx/gluon-${x}-100.txt \
#	${SAT}/Run2/runa-bjorx-4/gluon-${x}-100.txt 
	
#./plot.py -s "./kt-formula-report1/GBW-thresh-gluon-${x}-100.png" -y "linear" -c "b r r" -p "-. -. -" -l "GBW-r:GBW-kt:GBW-kt-run" \
#	${SAT}/Run2/GBW/gluon-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx2/gluon-${x}-100.txt \
#	${SAT}/Run2/runa-bjorx2-4/gluon-${x}-100.txt 
	
#./plot.py -s "./kt-formula-report1/BGK-gluon-${x}-100.png" -y "linear" -c "b r " -p "-. -" -l "BGK-r:BGK-kt" \
#	${SAT}/Run2/BGK/gluon-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx-BGK/gluon-${x}-100.txt 
	
	
#./plot.py -s "./kt-formula-report1/GBW-dipole-${x}-100.png" -c "b r r" -p "-. -. -" -l "GBW-r:GBW-kt:GBW-kt-run" \
#	${SAT}/Run2/GBW/dipole-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx/dipole-${x}-100.txt \
#	${SAT}/Run2/runa-bjorx-4/dipole-${x}-100.txt 
	
#./plot.py -s "./kt-formula-report1/GBW-thresh-dipole-${x}-100.png"  -c "b r r" -p "-. -. -" -l "GBW-r:GBW-kt:GBW-kt-run" \
#	${SAT}/Run2/GBW/dipole-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx2/dipole-${x}-100.txt \
#	${SAT}/Run2/runa-bjorx2-4/dipole-${x}-100.txt 
	
#./plot.py -s "./kt-formula-report1/BGK-dipole-${x}-100.png" -c "b r " -p "-. -" -l "BGK-r:BGK-kt" \
#	${SAT}/Run2/BGK/dipole-${x}-100.txt \
#	${SAT}/Run2/fixa-bjorx-BGK/dipole-${x}-100.txt 
#done
 
