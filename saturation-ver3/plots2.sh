#! /usr/bin/env bash
SAT=/home/tomoki/Saturation-Model/saturation-ver3

./plot.py -s "./kt-formula-report1/GBW-gluon.png" -y "linear" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	 -a "$ k_t^2\\; [\\mathrm{GeV^2}]$ :$\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx/gluon-2-100.txt \
	${SAT}/Run2/runa-bjorx-4/gluon-2-100.txt \
	${SAT}/Run2/GBW/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx/gluon-6-100.txt \
	${SAT}/Run2/runa-bjorx-4/gluon-6-100.txt \
	${SAT}/Run2/GBW/gluon-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-dipole.png" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	-a "$ r\\; [\\mathrm{GeV}]$ :$\\sigma/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx/dipole-2-100.txt \
	${SAT}/Run2/runa-bjorx-4/dipole-2-100.txt \
	${SAT}/Run2/GBW/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx/dipole-6-100.txt \
	${SAT}/Run2/runa-bjorx-4/dipole-6-100.txt \
	${SAT}/Run2/GBW/dipole-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-saturation.png" -a "x: $ Q_s^2(x)\\;[\\mathrm{GeV^2}]$" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	${SAT}/Run2/fixa-bjorx/saturation.txt \
	${SAT}/Run2/runa-bjorx-4/saturation.txt \
	${SAT}/Run2/GBW/saturation.txt 
	
./plot.py -s "./kt-formula-report1/GBW-thresh-gluon.png" -y "linear" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	-a "$ k_t^2\\; [\\mathrm{GeV^2}]$ :$\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx2/gluon-2-100.txt \
	${SAT}/Run2/runa-bjorx2-4/gluon-2-100.txt \
	${SAT}/Run2/GBW/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx2/gluon-6-100.txt \
	${SAT}/Run2/runa-bjorx2-4/gluon-6-100.txt \
	${SAT}/Run2/GBW/gluon-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-thresh-dipole.png" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	-a "$ r\\;[\\mathrm{GeV}]$ :$\\sigma/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx2/dipole-2-100.txt \
	${SAT}/Run2/runa-bjorx2-4/dipole-2-100.txt \
	${SAT}/Run2/GBW/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx2/dipole-6-100.txt \
	${SAT}/Run2/runa-bjorx2-4/dipole-6-100.txt \
	${SAT}/Run2/GBW/dipole-6-100.txt 
	
./plot.py -s "./kt-formula-report1/GBW-thresh-saturation.png" -a "x:$ Q_s^2(x)\\;[\\mathrm{GeV^2}]$" \
	-c "r r b r r b" -p "-. - -. -. - -." -l "GBW-kt:GBW-kt-run:GBW-r::::::" \
	${SAT}/Run2/fixa-bjorx2/saturation.txt \
	${SAT}/Run2/runa-bjorx2-4/saturation.txt\
	${SAT}/Run2/GBW/saturation.txt 

./plot.py -s "./kt-formula-report1/BGK-gluon.png" -y "linear" \
	-c "r b r b" -p "- -. - -." -l "BGK-kt:BGK-r::::::" \
	 -a "$ k_t^2 \\;[\\mathrm{GeV^2}] $ : $\\alpha\\mathcal{F}(x,k_t^2)/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx-BGK/gluon-2-100.txt \
	${SAT}/Run2/BGK/gluon-2-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/gluon-6-100.txt \
	${SAT}/Run2/BGK/gluon-6-100.txt 

./plot.py -s "./kt-formula-report1/BGK-dipole.png" \
	-c "r b r b" -p "- -. - -." -l "BGK-kt:BGK-r::::::" -a "r\\; [\\mathrm{GeV}]:$\\sigma/\\sigma_0$" \
	${SAT}/Run2/fixa-bjorx-BGK/dipole-2-100.txt \
	${SAT}/Run2/BGK/dipole-2-100.txt \
	${SAT}/Run2/fixa-bjorx-BGK/dipole-6-100.txt \
	${SAT}/Run2/BGK/dipole-6-100.txt 
	
	
./plot.py -s "./kt-formula-report1/BGK-saturation.png" \
	-c "r b r b" -p "- -. - -." -l "BGK-kt:BGK-r::::::" -a "x:$ Q_s^2(x)\\;[\\mathrm{GeV^2}]$" \
	${SAT}/Run2/fixa-bjorx-BGK/saturation.txt \
	${SAT}/Run2/BGK/saturation.txt 
	
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
 
