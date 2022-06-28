#!/usr/bin/env bash 

#mkdir PlotRelated

BGKSDIR=./BGKSStar/Mass0.0-Qup500-Model3-Sud1-rfix0
BGKDIR=./BGKStar/Mass0.0-Qup500-Model1-Sud0-rfix0
GBWDIR=./GBW/Mass0.0-Qup500-Model0-Sud0-rfix0
GBWSDIR=./GBWS/Mass0.0-Qup500-Model22-Sud1-rfix0

${BGKDIR}/dipole -in ${BGKSDIR}/result.txt -out PlotRelated/BGKdipole-500-5.txt -Q2 500 -x 5 
${GBWDIR}/dipole -in ${GBWSDIR}/result.txt -out PlotRelated/GBWdipole-500-5.txt -Q2 500 -x 5 

${BGKDIR}/F2-slope -in ${BGKSDIR}/result.txt -out PlotRelated/BGKF2-slope-3.txt -Q2 500 -x 3 
${GBWDIR}/F2-slope -in ${GBWSDIR}/result.txt -out PlotRelated/GBWF2-slope-3.txt -Q2 500 -x 3 


./plot.py \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ./PlotRelated/BGK-effect-of-sudakov.png \
 -a "r (fm):sigma/sigma_0"\
 ${BGKDIR}/dipole-500-5.txt \
 ${BGKSDIR}/dipole-500-5.txt \
 ./PlotRelated/BGKdipole-500-5.txt
 
./plot.py \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ./PlotRelated/GBW-effect-of-sudakov.png \
 -a "r (fm):sigma/sigma_0"\
 ${GBWDIR}/dipole-500-5.txt \
 ${GBWSDIR}/dipole-500-5.txt \
 ./PlotRelated/GBWdipole-500-5.txt

./plot.py \
 -y linear \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ./PlotRelated/BGK-effect-of-sudakov-slope.png \
 -a "Q^2 : slope"\
 ${BGKDIR}/F2-slope-3.txt \
 ${BGKSDIR}/F2-slope-3.txt \
 ./PlotRelated/BGKF2-slope-3.txt
 
./plot.py \
 -y linear \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ./PlotRelated/GBW-effect-of-sudakov-slope.png \
 -a "slope :Q^2"\
 ${GBWDIR}/F2-slope-3.txt \
 ${GBWSDIR}/F2-slope-3.txt \
 ./PlotRelated/GBWF2-slope-3.txt


 ./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "Results of Massless Fits" \
 -s ./PlotRelated/MasslessModels.png \
 -a "r (fm):sigma/sigma_0"\
 ${GBWDIR}/dipole-500-5.txt\
 ${GBWSDIR}/dipole-500-5.txt\
 ${BGKDIR}/dipole-500-5.txt \
 ${BGKSDIR}/dipole-500-5.txt 

 ./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "F2 -  Massless Fits" \
 -s ./PlotRelated/MasslessModelsF2.png \
 -a "Q^2 GeV:F2"\
 ${GBWDIR}/F2-3.txt\
 ${GBWSDIR}/F2-3.txt\
 ${BGKDIR}/F2-3.txt \
 ${BGKSDIR}/F2-3.txt

 ./plot.py \
 -y linear\
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "F2 slope -  Massless Fits" \
 -s ./PlotRelated/MasslessModelsF2-slope.png \
 -a "Q^2 GeV:- dlog(F2)/dlog(x)"\
 ${GBWDIR}/F2-slope-3.txt\
 ${GBWSDIR}/F2-slope-3.txt\
 ${BGKDIR}/F2-slope-3.txt \
 ${BGKSDIR}/F2-slope-3.txt 

./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "Critical Lines -  Massless Fits Q2=500" \
 -s ./PlotRelated/MasslessModelsCritical500.png \
 -a " x : Q^2_s "\
 ${GBWDIR}/critical-500.txt\
 ${GBWSDIR}/critical-500.txt\
 ${BGKDIR}/critical-500.txt \
 ${BGKSDIR}/critical-500.txt 


./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "Critical Lines -  Massless Fits Q2=50" \
 -s ./PlotRelated/MasslessModelsCritical50.png \
 -a " x : Q^2_s "\
 ${GBWDIR}/critical-50.txt\
 ${GBWSDIR}/critical-50.txt\
 ${BGKDIR}/critical-50.txt \
 ${BGKSDIR}/critical-50.txt 
