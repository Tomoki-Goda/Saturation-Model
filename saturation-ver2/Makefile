#########################################################
# before use do export DIR=<directory >			#
# like 							#
# for i in ./GBW/* ;  DIR=${i}; make ;done	#
#							#
#########################################################



DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h 

LIBS=-lm -lmathlib -lkernlib -lpacklib	 -Wno-unused-result

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR}



all: clean main critical dipole F2 F2-slope psi

fit: clean main

plot: critical dipole F2 F2-slope psi


main: ${DEPS} ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o 
	gcc ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/main -Ofast ${INC}

critical : ${DEPS} ${UTILDIR}/critical-line.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/critical-line.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/critical -Ofast  ${INC}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/dipole -Ofast  ${INC}

F2: ${DEPS} ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2 -Ofast  ${INC}

F2-slope: ${DEPS} ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2-slope -Ofast  ${INC}

psi: ${DEPS} ${UTILDIR}/Plot-Psi.c ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o
	gcc ${UTILDIR}/Plot-Psi.c  ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/psi -Ofast  ${INC}

curvature : ${DEPS} ${UTILDIR}/curvature.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o
	gcc ${UTILDIR}/curvature.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/curvature -Ofast  ${INC}



${DIR}/read-and-fit.o: ${DEPS} ${FUNCDIR}/read-and-fit.c 
	echo ${DIR}
	gcc -c ${FUNCDIR}/read-and-fit.c -Ofast  ${LIBS}  ${INC}  -o ${DIR}/read-and-fit.o

#sudakov.o: ${DEPS} sudakov.c mu2.h
#	gcc -c sudakov.c -Ofast  ${LIBS} -I ${DIR}  -o ${DIR}/sudakov.o

#gluon-chebyshev.o: ${DEPS} gluon-chebyshev.c gluons.h chebyshev-1.h
#	gcc -c gluon-chebyshev.c  ${LIBS} -Ofast -I ${DIR}  -o ${DIR}/gluon-chebyshev.o 
       
${DIR}/dipole-cross-section.o: ${DEPS} ${FUNCDIR}/dipole-cross-section.c
	gcc -c ${FUNCDIR}/dipole-cross-section.c  ${LIBS} -Ofast  ${INC}  -o ${DIR}/dipole-cross-section.o 

${DIR}/photon-wave-function.o: ${DEPS}  ${FUNCDIR}/photon-wave-function.c 
	gcc -c ${FUNCDIR}/photon-wave-function.c  ${LIBS} -Ofast  ${INC}  -o ${DIR}/photon-wave-function.o 
	
${DIR}/simpson-integral.o: ${FUNCDIR}/simpson-integral.c
	gcc -c ${FUNCDIR}/simpson-integral.c  ${LIBS} -Ofast   ${INC} -o ${DIR}/simpson-integral.o 
       	
clean: 
	rm -f ${DIR}/*.o
	rm -f ${DIR}/main
	 		
