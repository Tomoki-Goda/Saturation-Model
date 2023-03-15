#! /usr/bin/env bash
SAT=/home/tomoki/saturation-ver2/Run2
#SAT=/home/tomoki/Saturation-Model/saturation-ver3/Run2802


#./Plotting/surface.py -s "./kt-formula-report1/dipole-BGK.png"  -d \
#	-m ./Run3/BGK/dipole-grid.txt  ./Run3/fixa-bjorx-BGK/dipole-grid.txt #./Run3/fixa-bjorx-BGK2/dipole-grid.txt



./plot.py -s "./kt-formula-report1/Sudakov/GBWS-gluon.png" -y "linear" \
	-c "b r b r" -p "-. - -. -" -l "GBW:GBW+Sud::" \
	 -a "\$ k_t^2\\; [\\mathrm{GeV^2}]\$ :\$\\alpha\\mathcal{F}(x,k_t^2)\$" \
	${SAT}/GBW/Mass0.0-Qup650-Model0-Sud0/gluon-0-2.txt \
	${SAT}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/gluon-0-2.txt \
	${SAT}/GBW/Mass0.0-Qup650-Model0-Sud0/gluon-0-6.txt \
	${SAT}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/gluon-0-6.txt 
	
./plot.py -s "./kt-formula-report1/GBWS-dipole.png" \
	-c "b r  b r" -p "-. - -. -" -l "GBW:GBW+Sud::" \
	-a "\$ r\\; [\\mathrm{GeV}]\$ :\$\\sigma/\\sigma_0\$" \
	${SAT}/GBW/Mass0.0-Qup650-Model0-Sud0/dipole-0-2.txt \
	${SAT}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/dipole-0-2.txt \
	${SAT}/GBW/Mass0.0-Qup650-Model0-Sud0/dipole-0-6.txt \
	${SAT}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/dipole-0-6.txt 
	
./plot.py -s "./kt-formula-report1/Sudakov/GBWS-saturation.png" -a "x: \$ Q_s^2(x)\\;[\\mathrm{GeV^2}]\$" \
	-c "b r" -p "-. -" -l "GBW:GBW+Sud" \
	${SAT}/GBW/Mass0.0-Qup650-Model0-Sud0/tmd-critical-0.txt \
	${SAT}/GBWS-Fix-S/Mass0.0-Qup650-Model22-Sud1/tmd-critical-0.txt 
	
./plot.py -s "./kt-formula-report1/Sudakov/BGKS-gluon.png" -y "linear" \
	-c "b r b r" -p "-. - -. -" -l "BGK:BGK+Sud::" \
	 -a "\$ k_t^2\\; [\\mathrm{GeV^2}]\$ :\$\\alpha\\mathcal{F}(x,k_t^2)\$" \
	${SAT}/BGK/Mass0.0-Qup650-Model1-Sud0/gluon-0-2.txt \
	${SAT}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/gluon-0-2.txt \
	${SAT}/BGK/Mass0.0-Qup650-Model1-Sud0/gluon-0-6.txt \
	${SAT}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/gluon-0-6.txt 
	
./plot.py -s "./kt-formula-report1/Sudakov/BGKS-dipole.png" \
	-c "b r  b r" -p "-. - -. -" -l "BGK:BGK+Sud::" \
	-a "\$ r\\; [\\mathrm{GeV}]\$ :\$\\sigma/\\sigma_0\$" \
	${SAT}/BGK/Mass0.0-Qup650-Model1-Sud0/dipole-0-2.txt \
	${SAT}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/dipole-0-2.txt \
	${SAT}/BGK/Mass0.0-Qup650-Model1-Sud0/dipole-0-6.txt \
	${SAT}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/dipole-0-6.txt 
	
./plot.py -s "./kt-formula-report1/Sudakov/BGKS-saturation.png" -a "x: \$ Q_s^2(x)\\;[\\mathrm{GeV^2}]\$" \
	-c "b r" -p "-. -" -l "BGK:BGK+Sud" \
	${SAT}/BGK/Mass0.0-Qup650-Model1-Sud0/tmd-critical-0.txt \
	${SAT}/BGKS-Fix-S/Mass0.0-Qup650-Model3-Sud1/tmd-critical-0.txt 
	

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
 
