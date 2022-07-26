#! /usr/bin/env bash

RUNDIR=../Run2607

rm -r ${RUNDIR}/*/*

 ./Auto-Control -dir ${RUNDIR}-Fejer/GBW -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 0
 ./Auto-Control -dir ${RUNDIR}/GBW -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 0
 ./Auto-Control -dir ${RUNDIR}-Fejer/GBWS -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWS -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}-Fejer/GBWS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}-Fejer/GBWSnp-Fix-S -rfix 0 -sudakov 2 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWSnp-Fix-S -rfix 0 -sudakov 2 -lmass 0.0 0.0196 -qup 100 650 -model 22

 ./Auto-Control -dir ${RUNDIR}-Fejer/BGKS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 3
 ./Auto-Control -dir ${RUNDIR}/BGKS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 3
 ./Auto-Control -dir ${RUNDIR}-Fejer/BGK -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 1
 ./Auto-Control -dir ${RUNDIR}/BGK -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 1


 ./Append "#define FEJER 1" ${RUNDIR}-Fejer/*/M* 
 ./Append "#define INDEPENDENT_RMAX 1" ${RUNDIR}*/BGKS*/M*
# ./Append "#define N_CHEB_R 100" ${RUNDIR}-Fejer/*/M* 
# ./Append "#define N_SIMPS_R 50" ${RUNDIR}/*/M* 

#/////////////////////////////////////////////////////////////////////////////

export DIR

for i in ${RUNDIR}*/GBW/M*
do 
	cp ./minuit-strategy-gbw.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/GBWS-Fix-S/M*
do 
	cp ./minuit-strategy-gbws-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/GBW-Fix-C/M*
do 
	cp ./minuit-strategy-gbws-fix-c.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/GBWS-Fix-rmax/M*
do 
	cp ./minuit-strategy-gbws-fix-rmax.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/GBWS/M*
do 
	cp ./minuit-strategy-gbws.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/GBWSnp-Fix-S/M*
do 
	cp ./minuit-strategy-gbws-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make
done



for i in ${RUNDIR}*/BGK/M*
do 
	cp ./minuit-strategy-bgk.h ${i}/minuit-run.h
	DIR=${i}
	make
done

for i in ${RUNDIR}*/BGKS-Fix-S/M*
do 
	cp ./minuit-strategy-bgks-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make
done























