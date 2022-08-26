#! /usr/bin/env bash

RUNDIR=../Run2608

echo "rebuild ${RUNDIR} [y/n]?"
read ans
if [ ${ans} != 'y' ]
then
	echo ${ans}
	echo "abortig"	
	exit 1
fi

rm -r ${RUNDIR}/*

for i in GBW BGK GBWS-Fix-S BGKS-Fix-S GBWSnp-Fix-S BGKSnp-Fix-S; do mkdir ${RUNDIR}/${i} ; done 

 ./Auto-Control -dir ${RUNDIR}/GBWS-Fix-S -sudakov 1 -lmass 0.0 0.0196 -qup  5 25 50 100 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/BGKS-Fix-S -sudakov 1 -lmass 0.0 0.0196 -qup  5 25 50 100 650 -model 3
 ./Auto-Control -dir ${RUNDIR}/GBW -sudakov 0 -lmass 0.0 0.0196 -qup  5 25 50 100 650 -model 0
 ./Auto-Control -dir ${RUNDIR}/BGK -sudakov 0 -lmass 0.0 0.0196 -qup  5 25 50 100 650 -model 1

 ./Auto-Control -dir ${RUNDIR}/GBWSnp-Fix-S -sudakov 2 -lmass 0.0 0.0196 -qup 650 -model 22
 ./Auto-Control -dir ${RUNDIR}/BGKSnp-Fix-S -sudakov 2 -lmass 0.0 0.0196 -qup 650 -model 3


 ./Append "#define INDEPENDENT_RMAX 1" ${RUNDIR}*/BGKS*/M*

 for i in BGKS-Share-rmax BGKSnp-Share-rmax; do mkdir ${RUNDIR}/${i} ; done
 ./Auto-Control -dir ${RUNDIR}/BGKS-Share-rmax -sudakov 1 -lmass 0.0 0.0196 -qup  650 -model 22
 ./Auto-Control -dir ${RUNDIR}/BGKSnp-Share-rmax  -sudakov 2 -lmass 0.0 0.0196 -qup  650 -model 3
 
 
 ./Append "#define FEJER 1" ${RUNDIR}/*/M* 
 ./Append "#define N_CHEB_R 120" ${RUNDIR}/*/M* 
 ./Append "#define DGAUSS_PREC 1.0e-8" ${RUNDIR}/*/M* 
#/////////////////////////////////////////////////////////////////////////////

export DIR

for i in ${RUNDIR}*/GBW/M*
do 
	cp ./minuit-strategy-gbw.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/GBWS-Fix-S/M*
do 
	cp ./minuit-strategy-gbws-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/GBWSnp-Fix-S/M*
do 
	cp ./minuit-strategy-gbws-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done


#for i in ${RUNDIR}*/GBWS-Fix-C/M*
#do 
#	cp ./minuit-strategy-gbws-fix-c.h ${i}/minuit-run.h
#	DIR=${i}
#	make --silent
#done

######################################################################


for i in ${RUNDIR}*/BGK/M*
do 
	cp ./minuit-strategy-bgk.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/BGKS-Fix-S/M*
do 
	cp ./minuit-strategy-bgks-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/BGKSnp-Fix-S/M*
do 
	cp ./minuit-strategy-bgks-fix-s.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/BGKS-Share-rmax/M*
do 
	cp ./minuit-strategy-bgks-share-rmax.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done

for i in ${RUNDIR}*/BGKSnp-Share-rmax/M*
do 
	cp ./minuit-strategy-bgks-share-rmax.h ${i}/minuit-run.h
	DIR=${i}
	make --silent
done
#for i in ${RUNDIR}*/BGKS-Fix-C/M*
#do 
#	cp ./minuit-strategy-bgks-fix-c.h ${i}/minuit-run.h
#	DIR=${i}
#	make --silent
#done






















