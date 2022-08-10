#! /usr/bin/env bash

RUNDIR=../Run0308

rm -r ${RUNDIR}/*/*


for i in GBW \
       	GBWS \
	GBWS-Fix-S \
	GBWS-Fix-C \
	GBWS-Fix-rmax \
	GBWSnp-Fix-S \
	BGK \
	BGKS-Fix-S \
	BGKSnp-Fix-S		
do 
	mkdir ${RUNDIR}/${i}  
done


 ./Auto-Control -dir ${RUNDIR}/GBW -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 0
 ./Auto-Control -dir ${RUNDIR}/GBWS -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWSnp-Fix-S -rfix 0 -sudakov 2 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-C -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-rmax -rfix 0 -sudakov 2 -lmass 0.0 0.0196 -qup 100 650 -model 22

 ./Auto-Control -dir ${RUNDIR}/BGKSnp-Fix-S -rfix 0 -sudakov 2 -lmass 0.0 0.0196 -qup 100 650 -model 3
 ./Auto-Control -dir ${RUNDIR}/BGKS-Fix-S -rfix 0 -sudakov 1 -lmass 0.0 0.0196 -qup 100 650 -model 3
 ./Auto-Control -dir ${RUNDIR}/BGK -rfix 0 -sudakov 0 -lmass 0.0 0.0196 -qup 100 650 -model 1


 ./Append "#define FEJER 1" ${RUNDIR}/*/M* 
 ./Append "#define INDEPENDENT_RMAX 1" ${RUNDIR}*/BGKS*/M*
 ./Append "#define N_CHEB_R 80" ${RUNDIR}/*/M* 
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

for i in ${RUNDIR}*/GBWS-Fix-C/M*
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

for i in ${RUNDIR}*/BGKSnp-Fix-S/M*
do 
	cp ./minuit-strategy-bgksnp-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make
done























