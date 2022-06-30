#!/usr/bin/env bash 

#mkdir PlotRelated

BGKBGKSDIR=../BGK-BGKS/Mass0.0-Qup650-Model1-Sud0-rfix0
BGKSDIR=../BGKS/Mass0.0-Qup650-Model3-Sud1-rfix0
BGKDIR=../BGK/Mass0.0-Qup650-Model1-Sud0-rfix0
GBWGBWSDIR=../GBW-GBWS/Mass0.0-Qup650-Model0-Sud0-rfix0
GBWDIR=../GBW/Mass0.0-Qup650-Model0-Sud0-rfix0
GBWSDIR=./GBWS/Mass0.0-Qup650-Model22-Sud1-rfix0

#SAVEDIR=./PlotRelated
SAVEDIR="/media/tomoki/TOMOKI-USB/Saturation-Model/NewPlots"

#${BGKDIR}/dipole -in ${BGKSDIR}/result.txt -out PlotRelated/BGKdipole-500-5.txt -Q2 500 -x 5 
#${GBWDIR}/dipole -in ${GBWSDIR}/result.txt -out PlotRelated/GBWdipole-500-5.txt -Q2 500 -x 5 

#${BGKDIR}/F2-slope -in ${BGKSDIR}/result.txt -out PlotRelated/BGKF2-slope-3.txt -Q2 500 -x 3 
#${GBWDIR}/F2-slope -in ${GBWSDIR}/result.txt -out PlotRelated/GBWF2-slope-3.txt -Q2 500 -x 3 

############### Sudakov ###########

./plot.py \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov.png \
 -a "r (fm):sigma/sigma_0"\
 ${BGKDIR}/dipole-500-5.txt \
 ${BGKSDIR}/dipole-500-5.txt \
 ${BGKBGKSDIR}/dipole-500-5.txt 
 
./plot.py \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov.png \
 -a "r (fm):sigma/sigma_0"\
 ${GBWDIR}/dipole-500-5.txt \
 ${GBWSDIR}/dipole-500-5.txt \
 ${GBWGBWSDIR}/dipole-500-5.txt 

./plot.py \
 -y linear \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov-linear.png \
 -a "r (fm):sigma/sigma_0"\
 ${BGKDIR}/dipole-500-5.txt \
 ${BGKSDIR}/dipole-500-5.txt \
 ${BGKBGKSDIR}/dipole-500-5.txt 
 
./plot.py \
 -y linear \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov-linear.png \
 -a "r (fm):sigma/sigma_0"\
 ${GBWDIR}/dipole-500-5.txt \
 ${GBWSDIR}/dipole-500-5.txt \
 ${GBWGBWSDIR}/dipole-500-5.txt 

./plot.py \
./plot.py \
 -y linear \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov-slope.png \
 -a "Q^2 : slope"\
 ${BGKDIR}/F2-slope-3.txt \
 ${BGKSDIR}/F2-slope-3.txt \
 ${BGKBGKSDIR}/F2-slope-3.txt 
 
./plot.py \
 -y linear \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov-slope.png \
 -a "Q^2 : slope"\
 ${GBWDIR}/F2-slope-3.txt \
 ${GBWSDIR}/F2-slope-3.txt \
 ${GBWGBWSDIR}/F2-slope-3.txt 

./plot.py \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov-critical.png \
 -a "x : saturation scale"\
 ${BGKDIR}/critical-500.txt \
 ${BGKSDIR}/critical-500.txt \
 ${BGKBGKSDIR}/critical-500.txt 
 
./plot.py \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov-critical.png \
 -a "x : saturation scale"\
 ${GBWDIR}/critical-500.txt \
 ${GBWSDIR}/critical-500.txt \
 ${GBWGBWSDIR}/critical-500.txt 

./plot.py \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov-F2-3.png \
 -a "Q^2:F2"\
 ${GBWDIR}/F2-3.txt \
 ${GBWSDIR}/F2-3.txt \
 ${GBWGBWSDIR}/F2-3.txt 

./plot.py \
 -p "- -. --" \
 -l "GBW :GBWS :GBW with GBWS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/GBW-effect-of-sudakov-F2-5.png \
 -a "Q^2:F2"\
 ${GBWDIR}/F2-5.txt \
 ${GBWSDIR}/F2-5.txt \
 ${GBWGBWSDIR}/F2-5.txt 

./plot.py \
 -p "- -. --" \
 -l "BGK :BGKS :GBK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov-F2-3.png \
 -a "Q^2:F2"\
 ${BGKDIR}/F2-3.txt \
 ${BGKSDIR}/F2-3.txt \
 ${BGKBGKSDIR}/F2-3.txt 

./plot.py \
 -p "- -. --" \
 -l "BGK :BGKS :BGK with BGKS prameters" \
 -c "cyan red blue" \
 -t "Effect of Sudakov Factor" \
 -s ${SAVEDIR}/BGK-effect-of-sudakov-F2-5.png \
 -a "Q^2: F2"\
 ${BGKDIR}/F2-5.txt \
 ${BGKSDIR}/F2-5.txt \
 ${BGKBGKSDIR}/F2-5.txt 



###########   Models  ################
 ./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "Results of Massless Fits" \
 -s ${SAVEDIR}/MasslessModels.png \
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
 -s ${SAVEDIR}/MasslessModelsF2.png \
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
 -s ${SAVEDIR}/MasslessModelsF2-slope.png \
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
 -s ${SAVEDIR}/MasslessModelsCritical500.png \
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
 -s ${SAVEDIR}/MasslessModelsCritical50.png \
 -a " x : Q^2_s "\
 ${GBWDIR}/critical-50.txt\
 ${GBWSDIR}/critical-50.txt\
 ${BGKDIR}/critical-50.txt \
 ${BGKSDIR}/critical-50.txt 
 
 
./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "F2 -  Massless Fits Q2=500 x=10^-3" \
 -s ${SAVEDIR}/MasslessModelsF2-3.png \
 -a "Q^2 :F2 "\
 ${GBWDIR}/F2-3.txt\
 ${GBWSDIR}/F2-3.txt\
 ${BGKDIR}/F2-3.txt \
 ${BGKSDIR}/F2-3.txt 
 
./plot.py \
 -p "-- -. -- -." \
 -l "GBW :GBWS :BGK :BGKS" \
 -c "green orange red blue" \
 -t "F2 -  Massless Fits Q2=500 x=10^-5" \
 -s ${SAVEDIR}/MasslessModelsF2-5.png \
 -a "Q^2 :F2 "\
 ${GBWDIR}/F2-5.txt\
 ${GBWSDIR}/F2-5.txt\
 ${BGKDIR}/F2-5.txt \
 ${BGKSDIR}/F2-5.txt 
 
 
#################### Q dep 
  ./plot.py \
  -y linear\
  -p "-- -. -- -. " \
  -l "GBW 50 GeV:GBW 500 GeV:GBWS 50 GeV:GBWS 500 GeV" \
  -c "cyan magenta blue red" \
  -t "Q dependence of dipole cross-section -  Massless Fits Q_max=650" \
  -s ${SAVEDIR}/GBWQdependence.png \
  -a " r fm : sigma/sigma_0 "\
   ./GBW/Mass0.0-Qup650-Model0-Sud0-rfix0/dipole-50-5.txt \
   ./GBW/Mass0.0-Qup650-Model0-Sud0-rfix0/dipole-500-5.txt \
   ./GBWS/Mass0.0-Qup650-Model22-Sud1-rfix0/dipole-50-5.txt \
   ./GBWS/Mass0.0-Qup650-Model22-Sud1-rfix0/dipole-500-5.txt 
   
    ./plot.py \
  -y linear\
   -p "-- -. -- -." \
  -l "BGK 50 GeV:BGK 500 GeV:BGKS 50 GeV:BGKS 500 GeV" \
  -c "cyan magenta blue red" \
  -t "Q dependence of dipole cross-section -  Massless Fits Q_max=650" \
  -s ${SAVEDIR}/BGKQdependence.png \
  -a " r fm : sigma/sigma_0 "\
   ./BGK/Mass0.0-Qup650-Model1-Sud0-rfix0/dipole-50-5.txt \
   ./BGK/Mass0.0-Qup650-Model1-Sud0-rfix0/dipole-500-5.txt \
   ./BGKS/Mass0.0-Qup650-Model3-Sud1-rfix0/dipole-50-5.txt \
   ./BGKS/Mass0.0-Qup650-Model3-Sud1-rfix0/dipole-500-5.txt 

