#########################################################
# before use do export DIR=<directory >			#
# like 							#
# for i in ./GBW/* ;  DIR=${i}; make ;done	#
#							#
#########################################################



DEPS= ${DIR}/control.h  ./control-default.h ./constants.h 


LIBS=-lm -lmathlib -lkernlib -lpacklib	 -Wno-unused-result

UTILDIR=./Utilities


all: main critical dipole F2 F2-slope psi
	echo ${DIR}



main: ${DEPS} main.c read-and-fit.o sudakov.o gluon-chebyshev.o dipole-cross-section.o photon-wave-function.o simpson-integral.o 
	gcc main.c ${DIR}/read-and-fit.o ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/main -Ofast -I ${DIR}

critical : ${DEPS} ${UTILDIR}/critical-line.c sudakov.o gluon-chebyshev.o dipole-cross-section.o
	gcc ${UTILDIR}/critical-line.c  ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/critical -Ofast -I ${DIR}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c sudakov.o gluon-chebyshev.o dipole-cross-section.o
	gcc ${UTILDIR}/Plot-DP.c ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/dipole -Ofast -I ${DIR}

F2: ${DEPS} ${UTILDIR}/Plot-F2.c sudakov.o gluon-chebyshev.o dipole-cross-section.o photon-wave-function.c simpson-integral.c
	gcc ${UTILDIR}/Plot-F2.c ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2 -Ofast -I ${DIR}

F2-slope: ${DEPS} ${UTILDIR}/F2-slope.c sudakov.o gluon-chebyshev.o dipole-cross-section.o photon-wave-function.c simpson-integral.c
	gcc ${UTILDIR}/F2-slope.c ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2-slope -Ofast -I ${DIR}

psi: ${DEPS} ${UTILDIR}/Plot-Psi.c photon-wave-function.c simpson-integral.c
	gcc ${UTILDIR}/Plot-Psi.c  ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/psi -Ofast -I ${DIR}



read-and-fit.o: ${DEPS} read-and-fit.c 
	gcc -c read-and-fit.c -Ofast  ${LIBS} -I ${DIR}  -o ${DIR}/read-and-fit.o

sudakov.o: ${DEPS} sudakov.c mu2.h
	gcc -c sudakov.c -Ofast  ${LIBS} -I ${DIR}  -o ${DIR}/sudakov.o

gluon-chebyshev.o: ${DEPS} gluon-chebyshev.c gluons.h chebyshev-1.h
	gcc -c gluon-chebyshev.c  ${LIBS} -Ofast -I ${DIR}  -o ${DIR}/gluon-chebyshev.o 
       
dipole-cross-section.o: ${DEPS} dipole-cross-section.c
	gcc -c dipole-cross-section.c  ${LIBS} -Ofast -I ${DIR}  -o ${DIR}/dipole-cross-section.o 

photon-wave-function.o: ${DEPS}  photon-wave-function.c 
	gcc -c photon-wave-function.c  ${LIBS} -Ofast -I ${DIR}  -o ${DIR}/photon-wave-function.o 
	
simpson-integral.o: simpson-integral.c
	gcc -c simpson-integral.c  ${LIBS} -Ofast   -o ${DIR}/simpson-integral.o 
       	
clean: 
	rm ${DIR}/*.o
	rm ${DIR}/main
	 		
