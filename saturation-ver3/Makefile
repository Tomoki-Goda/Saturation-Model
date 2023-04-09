#########################################################
# before use do export DIR=<directory >			#
# like 							#
# export DIR						#
# for i in ./GBW/* ;  DIR=${i}; make ;done	        #
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.hh 
LIBS=-lm -lcuba -lMinuit2 -lgsl -lgslcblas  -fopenmp

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR} 

#OPT= -Wno-unused-result -fp-model=precise -O3#-g # -Wno-format 
OPT=-O3 -Wno-unused-result 

#CXX=/opt/intel/oneapi/compiler/2023.0.0/linux/bin/icpx
#CXX=icpx
CXX=g++

########################################################################################################
all: main

main: ${DIR}/gluons.o ${DIR}/r-formula.o ${DIR}/gluon-integrand.o   ${DIR}/dipole-gluon.o \
		${DIR}/interpolation-gluon.o ${DIR}/main.o
	${CXX} $^ -o ${DIR}/main ${INC} ${OPT} ${LIBS} -I${FUNCDIR} 

${DIR}/%.o : ${FUNCDIR}/%.cc 
		echo $@
		${CXX} -c -I${INCDIR} -I${FUNCDIR} -I${DIR} ${OPT} ${LIBS} $< -o $@
#.cc.o:
#	${CXX} -c -I${INCDIR} -I${FUNCDIR}  ${FUNCDIR}/$<
clean: 
	rm -f ${DIR}/*.o
	rm -f ${FUNCDIR}/*.o
	rm -f ${DIR}/main
	 		
