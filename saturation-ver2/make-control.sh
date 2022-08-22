#! /usr/bin/env bash

RUNDIR=../Run1708

rm -r ${RUNDIR}/*/*


for i in GBW GBWS-Fix-S BGK BGKS-Fix-S; do mkdir ${RUNDIR}/${i} ; done 

 ./Auto-Control -dir ${RUNDIR}/GBW -sudakov 0 -lmass 0.0 0.0196 -qup 5 10 25 50 75 100 650 -model 0
 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-S -sudakov 1 2 -lmass 0.0 0.0196 -qup 5 10 25 50 75 100 650 -model 22

 ./Auto-Control -dir ${RUNDIR}/BGK -sudakov 0 -lmass 0.0 0.0196 -qup 5 10 25 50 75 100 650 -model 1
 ./Auto-Control -dir ${RUNDIR}/BGKS-Fix-S -sudakov 1 2 -lmass 0.0 0.0196 -qup 5 10 25 50 75 100 650 -model 3


 ./Append "#define FEJER 1" ${RUNDIR}/*/M* 
 ./Append "#define INDEPENDENT_RMAX 1" ${RUNDIR}*/BGKS*/M*
 ./Append "#define N_CHEB_R 120" ${RUNDIR}/*/M* 
 ./Append "#define DGAUSS_PREC 1.0e-5" ${RUNDIR}/*/M* 
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

#for i in ${RUNDIR}*/GBWS-Fix-C/M*
#do 
#	cp ./minuit-strategy-gbws-fix-c.h ${i}/minuit-run.h
#	DIR=${i}
#	make
#done

#for i in ${RUNDIR}*/GBWS-Fix-rmax/M*
#do 
#	cp ./minuit-strategy-gbws-fix-rmax.h ${i}/minuit-run.h
#	DIR=${i}
#	make
#done

#for i in ${RUNDIR}*/GBWS/M*
#do 
#	cp ./minuit-strategy-gbws.h ${i}/minuit-run.h
#	DIR=${i}
#	make
#done

#for i in ${RUNDIR}*/GBWSnp-Fix-S/M*
#do 
#	cp ./minuit-strategy-gbws-fix-s.h ${i}/minuit-run.h
#	DIR=${i}
#	make
#done



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























