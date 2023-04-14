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
all: clean  main

main: $(DIR)/gluons.o ${DIR}/r-formula.o ${DIR}/gluon-integrand.o   ${DIR}/dipole-gluon.o \
		${DIR}/interpolation-gluon.o ${DIR}/main.o 
	${CXX} $^ -o ${DIR}/main ${INC} ${OPT} ${LIBS} -I${FUNCDIR} 


$(DIR)/%.o: ${FUNCDIR}/%.cc ${FUNCDIR}/%.hh ${DEPS} 
		${CXX} -c -I${INCDIR} -I${FUNCDIR} -I${DIR} ${OPT} ${LIBS} $< -o $@


gluonsf.o:${FUNCDIR}/gluons.f
		gfortran -c -I${INCDIR} -I${FUNCDIR} -I${FUNCDIR}/fortran ${OPT} ${LIBS} $< -o $@
gluonsc.o:${FUNCDIR}/gluons.cc
		${CXX} -c -I${INCDIR} -I${FUNCDIR}  ${OPT} ${LIBS} $< -o $@

test: gluonsf.o gluonsc.o
	${CXX} -c test-gluons.cc -I${FUNCDIR}  -o test-gluons.o
	${CXX} test-gluons.o -o test-gluons -I${INCDIR} -I${FUNCDIR} -I. -I${DIR} ${OPT} ${LIBS} $^ -lgfortran -lmathlib

clean: 
	rm -f ${DIR}/*.o
	rm -f ${FUNCDIR}/*.o
	rm -f ${DIR}/main
	 		
