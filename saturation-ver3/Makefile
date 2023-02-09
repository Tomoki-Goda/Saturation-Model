#########################################################
# before use do export DIR=<directory >			#
# like 							#
# export DIR						#
# for i in ./GBW/* ;  DIR=${i}; make ;done	        #
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.hh 

LIBS=-lm -lmathlib -lkernlib -lpacklib -lquadmath -lcuba -lMinuit2 -lgsl -lgslcblas -lgfortran

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR} 

OPT=-O3 -Wno-unused-result -std=c++17#-g # -Wno-format 

#SUDDEP=./Functions/mu2.h ./Functions/sudakov.h
#GLUDEP=./Functions/gluon-chebyshev.h ./Functions/gluons.h ./Functions/chebyshev-1.h ./Functions/complex.h
#PARNIDIR=/home/tomoki/Tools/Parni

########################################################################################################
all: clean main 
grid: gluon-grid

fit: clean main

main: ${DEPS} ${FUNCDIR}/main.cc ${FUNCDIR}/kt-formula.hh ${FUNCDIR}/fcn.h 
	g++ ${FUNCDIR}/main.cc -o ${DIR}/main   ${INC} ${OPT} ${LIBS} 

result: ${DEPS} ${FUNCDIR}/format-result.cc ${FUNCDIR}/kt-formula.hh ${FUNCDIR}/fcn.h 
	g++ ${FUNCDIR}/format-result.cc ${LIBS} -o ${DIR}/result  ${INC} ${OPT}

#gluons.o: ${FUNCDIR}/gluons.f ${FUNCDIR}/clenshaw.f
#	gfortran ${FUNCDIR}/gluons.f -c -o ${FUNCDIR}/gluons.o

#main: ${DEPS} ${FUNCDIR}/main.cc ${DIR}/kt-formula.o 
#	g++ ${FUNCDIR}/main.cc  ${DIR}/kt-formula.o ${LIBS} -o ${DIR}/main  ${INC} ${OPT}

#${DIR}/kt-formula.o: ${FUNCDIR}/kt-formula.cc
#	g++ ${FUNCDIR}/kt-formula.cc -o ${DIR}/kt-formula.o -c  ${INC} ${OPT} 

plot: ${FUNCDIR}/plot-things.cc #${FUNCDIR}/gluons.o
	g++ ${FUNCDIR}/plot-things.cc -o ${DIR}/plot ${OPT} ${INC} ${LIBS}
	g++ ${FUNCDIR}/plot-grid.cc -o ${DIR}/grid  ${OPT} ${INC} ${LIBS}
#plot: dipole gluon

gluon: ${DEPS} ${UTILDIR}/tmd-gluon.c ${UTILDIR}/tmd-gluon-2.h
	g++ ${UTILDIR}/tmd-gluon.c -o ${DIR}/gluon ${OPT} ${INC} ${LIBS}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c 
	g++ ${UTILDIR}/Plot-DP.c -o ${DIR}/dipole ${OPT} ${INC} ${LIBS}

clean: 
	rm -f ${DIR}/*.o
	rm -f ${FUNCDIR}/*.o
	rm -f ${DIR}/main
	 		
