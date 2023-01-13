#########################################################
# before use do export DIR=<directory >			#
# like 							#
# export DIR						#
# for i in ./GBW/* ;  DIR=${i}; make ;done	        #
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.hh 

LIBS=-lm -lmathlib -lkernlib -lpacklib -lquadmath -lcuba -lMinuit2

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR} 

OPT=-O3 -Wno-unused-result #-g # -Wno-format 

SUDDEP=./Functions/mu2.h ./Functions/sudakov.h
GLUDEP=./Functions/gluon-chebyshev.h ./Functions/gluons.h ./Functions/chebyshev-1.h ./Functions/complex.h
PARNIDIR=/home/tomoki/Tools/Parni

########################################################################################################
all: clean main 
grid: gluon-grid

fit: clean main

main: ${DEPS} ${FUNCDIR}/main.cc ${DIR}/kt-formula.o ${PARNI}
	g++ ${FUNCDIR}/main.cc  ${DIR}/kt-formula.o ${PARNI} ${LIBS} -o ${DIR}/main  ${INC} ${OPT}

${DIR}/kt-formula.o: ${FUNCDIR}/kt-formula.cc
	g++ ${FUNCDIR}/kt-formula.cc -o ${DIR}/kt-formula.o -c  ${INC} ${OPT} 


		
clean: 
	rm -f ${DIR}/*.o
	rm -f ${DIR}/main
	 		
