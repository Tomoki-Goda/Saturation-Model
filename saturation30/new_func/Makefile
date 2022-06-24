#########################################################
# before use do export DIR=<directory >			#
# like 							#
# for i in ./GBW/* ;  DIR=${i}; make ;done	#
#							#
#########################################################



DEPS= ${DIR}/control.h  ./control-default.h ./constants.h 

full: main clean

all: setdir main
	echo ${DIR}


LIBS=-lm -lmathlib -lkernlib -lpacklib	 


main: ${DEPS} main.c read-and-fit.o sudakov.o gluon-chebyshev.o dipole-cross-section.o photon-wave-function.o simpson-integral.o 
	gcc main.c ${DIR}/read-and-fit.o ${DIR}/sudakov.o ${DIR}/gluon-chebyshev.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/main -Ofast -I ${DIR}

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


	 		
