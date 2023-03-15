#! /usr/bin/env bash
SAT=/home/tomoki/Saturation-Model/saturation-ver3

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW.png" \
	 -m ./Run3/GBW/gluon-grid.txt ./Run3/fixa-bjorx/gluon-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run.png"  \
	-m ./Run3/GBW/gluon-grid.txt ./Run3/runa-bjorx-4/gluon-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK.png"  \
	-m ./Run3/BGK/gluon-grid.txt ./Run3/fixa-bjorx-BGK/gluon-grid.txt #./Run3/fixa-bjorx-BGK2/gluon-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-contour.png"  \
	-m ./Run3/GBW/gluon-grid.txt ./Run3/fixa-bjorx/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-GBW-run-contour.png" \
	-m ./Run3/GBW/gluon-grid.txt ./Run3/runa-bjorx-4/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/gluon-BGK-contour.png" -c \
	-m ./Run3/BGK/gluon-grid.txt ./Run3/fixa-bjorx-BGK/gluon-grid.txt #./Run3/fixa-bjorx-BGK2/gluon-grid.txt -c

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW.png"  -d \
	-m ./Run3/GBW/dipole-grid.txt ./Run3/fixa-bjorx/dipole-grid.txt 

./Plotting/surface.py -s "./kt-formula-report1/dipole-GBW-run.png" -d \
	-m ./Run3/GBW/dipole-grid.txt ./Run3/runa-bjorx-4/dipole-grid.txt

./Plotting/surface.py -s "./kt-formula-report1/dipole-BGK.png"  -d \
	-m ./Run3/BGK/dipole-grid.txt  ./Run3/fixa-bjorx-BGK/dipole-grid.txt #./Run3/fixa-bjorx-BGK2/dipole-grid.txt



./plot.py -s "./kt-formula-report1/GBW-gluon.png" -y "linear" \
	-c "r r b r r b" -p "-. - -. -. - -." -l " \$k_t\$ -GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
	 -a "\$ k_t^2\\; [\\mathrm{GeV^2}]\$ :\$\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0\$" \
	${SAT}/Run3/fixa-bjorx/gluon-2-100.txt \
	${SAT}/Run3/runa-bjorx-4/gluon-2-100.txt \
	${SAT}/Run3/GBW/gluon-2-100.txt \
	${SAT}/Run3/fixa-bjorx/gluon-6-100.txt \
	${SAT}/Run3/runa-bjorx-4/gluon-6-100.txt \
	${SAT}/Run3/GBW/gluon-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-dipole.png" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "\$k_t\$-GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
	-a "\$ r\\; [\\mathrm{GeV}]\$ :\$\\sigma/\\sigma_0\$" \
	${SAT}/Run3/fixa-bjorx/dipole-2-100.txt \
	${SAT}/Run3/runa-bjorx-4/dipole-2-100.txt \
	${SAT}/Run3/GBW/dipole-2-100.txt \
	${SAT}/Run3/fixa-bjorx/dipole-6-100.txt \
	${SAT}/Run3/runa-bjorx-4/dipole-6-100.txt \
	${SAT}/Run3/GBW/dipole-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-saturation.png" -a "x: \$ Q_s^2(x)\\;[\\mathrm{GeV^2}]\$" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "\$k_t\$-GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
	${SAT}/Run3/fixa-bjorx/saturation.txt \
	${SAT}/Run3/runa-bjorx-4/saturation.txt \
	${SAT}/Run3/GBW/saturation.txt 
	
#./plot.py -s "./kt-formula-report1/GBW-thresh-gluon.png" -y "linear" \
#	-c "r r b r r b" -p "-. - -. -. - -." -l "\$k_t\$-GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
#	-a "\$ k_t^2\\; [\\mathrm{GeV^2}]\$ :\$\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0\$" \
#	${SAT}/Run3/fixa-bjorx2/gluon-2-100.txt \
#	${SAT}/Run3/runa-bjorx2-4/gluon-2-100.txt \
#	${SAT}/Run3/GBW/gluon-2-100.txt \
#	${SAT}/Run3/fixa-bjorx2/gluon-6-100.txt \
#	${SAT}/Run3/runa-bjorx2-4/gluon-6-100.txt \
#	${SAT}/Run3/GBW/gluon-6-100.txt 
	
#./plot.py -s "./kt-formula-report1/GBW-thresh-dipole.png" \
#	-c "r r b r r b" -p "-. - -. -. - -." -l "\$k_t\$-GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
#	-a "\$ r\\;[\\mathrm{GeV}]\$ :\$\\sigma/\\sigma_0\$" \
#	${SAT}/Run3/fixa-bjorx2/dipole-2-100.txt \
#	${SAT}/Run3/runa-bjorx2-4/dipole-2-100.txt \
#	${SAT}/Run3/GBW/dipole-2-100.txt \
#	${SAT}/Run3/fixa-bjorx2/dipole-6-100.txt \
#	${SAT}/Run3/runa-bjorx2-4/dipole-6-100.txt \
#	${SAT}/Run3/GBW/dipole-6-100.txt 
	
#./plot.py -s "./kt-formula-report1/GBW-thresh-saturation.png" -a "x:\$ Q_s^2(x)\\;[\\mathrm{GeV^2}]\$" \
#	-c "r r b r r b" -p "-. - -. -. - -." -l "\$k_t\$-GBW:rc-\$k_t\$-GBW:\$r\$-GBW::::::" \
#	${SAT}/Run3/fixa-bjorx2/saturation.txt \
#	${SAT}/Run3/runa-bjorx2-4/saturation.txt\
#	${SAT}/Run3/GBW/saturation.txt 

./plot.py -s "./kt-formula-report1/BGK-gluon.png" -y "linear" \
	-c "r b r b" -p "- -. - -." -l "\$k_t\$-BGK:\$r\$-BGK::::::" \
	 -a "\$ k_t^2 \\;[\\mathrm{GeV^2}] \$ : \$\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0\$" \
	${SAT}/Run3/fixa-bjorx-BGK/gluon-2-100.txt \
	${SAT}/Run3/BGK/gluon-2-100.txt \
	${SAT}/Run3/fixa-bjorx-BGK/gluon-6-100.txt \
	${SAT}/Run3/BGK/gluon-6-100.txt 

./plot.py -s "./kt-formula-report1/BGK-dipole.png" \
	-a "r\\; [\\mathrm{GeV}]:\$\\sigma/\\sigma_0\$" \
	-c "r b r b" -p "- -. - -." -l "\$k_t\$-BGK:\$r\$-BGK::::::" \
	${SAT}/Run3/fixa-bjorx-BGK/dipole-2-100.txt \
	${SAT}/Run3/BGK/dipole-2-100.txt \
	${SAT}/Run3/fixa-bjorx-BGK/dipole-6-100.txt \
	${SAT}/Run3/BGK/dipole-6-100.txt 
	
	
./plot.py -s "./kt-formula-report1/BGK-saturation.png" \
	-a "x:\$ Q_s^2(x)\\;[\\mathrm{GeV^2}]\$" \
	-c "r b r b" -p "- -. - -." -l "\$k_t\$-BGK:\$r\$-BGK::::::" \
	${SAT}/Run3/fixa-bjorx-BGK/saturation.txt \
	${SAT}/Run3/BGK/saturation.txt 
	
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
 
