#########################################################
# before use do export DIR=<directory >			#
# like 							#
# for i in ./GBW/* ;  DIR=${i}; make ;done	#
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.h 

LIBS=-lm -lmathlib -lkernlib -lpacklib  -Wno-unused-result

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR}

OPT=-O3

SUDDEP=./Functions/mu2.h ./Functions/sudakov.h
GLUDEP=./Functions/gluon-chebyshev.h ./Functions/gluons.h ./Functions/chebyshev-1.h ./Functions/complex.h

########################################################################################################

all: clean main critical dipole F2 F2-slope psi

fit: clean main

plot: critical dipole F2 F2-slope psi


main: ${DEPS} ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o 
	gcc ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/main ${INC} ${OPT} 



critical : ${DEPS} ${UTILDIR}/critical-line.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/critical-line.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/critical   ${INC} ${OPT}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/dipole  ${INC} ${OPT}

F2: ${DEPS} ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2   ${INC} ${OPT}

F2-slope: ${DEPS} ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2-slope   ${INC} ${OPT}

psi: ${DEPS} ${UTILDIR}/Plot-Psi.c ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/Plot-Psi.c  ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/psi   ${INC} ${OPT}

curvature : ${DEPS} ${UTILDIR}/curvature.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/curvature.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/curvature  ${INC} ${OPT}




${DIR}/read-and-fit.o: ${DEPS} ${FUNCDIR}/read-and-fit.c 
	echo ${DIR}
	gcc -c ${FUNCDIR}/read-and-fit.c   ${LIBS}  ${INC}  -o ${DIR}/read-and-fit.o ${OPT}

       
${DIR}/dipole-cross-section.o: ${DEPS} ${FUNCDIR}/dipole-cross-section.c ${GLUDEP} ${SUDDEP}
	gcc -c ${FUNCDIR}/dipole-cross-section.c  ${LIBS}   ${INC}  -o ${DIR}/dipole-cross-section.o  ${OPT}


${DIR}/photon-wave-function.o: ${DEPS}  ${FUNCDIR}/photon-wave-function.c 
	gcc -c ${FUNCDIR}/photon-wave-function.c  ${LIBS}   ${INC}  -o ${DIR}/photon-wave-function.o  ${OPT}


${DIR}/simpson-integral.o: ${FUNCDIR}/simpson-integral.c
	gcc -c ${FUNCDIR}/simpson-integral.c  ${LIBS}   ${INC} -o ${DIR}/simpson-integral.o  ${OPT}

		
clean: 
	rm -f ${DIR}/*.o
	rm -f ${DIR}/main
	 		
