#########################################################
# before use do export DIR=<directory >			#
# like 							#
# export DIR						#
# for i in ./GBW/* ;  DIR=${i}; make ;done	        #
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.h 

LIBS=-lm -lmathlib -lkernlib -lpacklib -lquadmath

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR}

OPT=-O3 -Wno-unused-result 

SUDDEP=./Functions/mu2.h ./Functions/sudakov.h
GLUDEP=./Functions/gluon-chebyshev.h ./Functions/gluons.h ./Functions/chebyshev-1.h ./Functions/complex.h

########################################################################################################

all: clean main critical dipole F2 FL F2-slope tmd-gluon tmd-gluon-x tmd-critical fcn export-data ww-gluon 
grid: gluon-grid

fit: clean main

plot: critical dipole F2 FL F2-slope tmd-gluon tmd-gluon-x tmd-critical fcn export-data  ww-gluon gluon-grid

main: ${DEPS} ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${FUNCDIR}/clenshaw.h 
	gcc ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/main ${INC} ${OPT}

gluon-grid: ${DEPS} ${UTILDIR}/gluon-grid.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon-2.h  
	gcc ${UTILDIR}/gluon-grid.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/gluon-grid   ${INC} ${OPT}

dipole-grid: ${DEPS} ${UTILDIR}/dipole-grid.c ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/dipole-grid.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/dipole-grid  ${INC} ${OPT}

critical : ${DEPS} ${UTILDIR}/critical-line.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/critical-line.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/critical   ${INC} ${OPT}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/dipole  ${INC} ${OPT}

F2: ${DEPS} ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o  ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/F2   ${INC} ${OPT}

FL: ${DEPS} ${UTILDIR}/Plot-FL.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/Plot-FL.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/FL   ${INC} ${OPT}

F2-slope: ${DEPS} ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/F2-slope   ${INC} ${OPT}

curvature : ${DEPS} ${UTILDIR}/curvature.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/curvature.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/curvature  ${INC} ${OPT}

tmd-gluon: ${DEPS} ${UTILDIR}/tmd-gluon.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon-2.h ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/tmd-gluon.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/tmd-gluon ${INC} ${OPT}

tmd-gluon-x: ${DEPS} ${UTILDIR}/tmd-gluon-x.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon-2.h ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/tmd-gluon-x.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/tmd-gluon-x   ${INC} ${OPT}

tmd-critical: ${DEPS} ${UTILDIR}/tmd-critical.c  ${UTILDIR}/tmd-gluon-2.h ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/tmd-critical.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/tmd-critical   ${INC} ${OPT}

ww-gluon: ${DEPS} ${UTILDIR}/ww-gluon.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/ww-gluon.c ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/ww-gluon   ${INC} ${OPT}

fcn: ${DEPS} ${FUNCDIR}/fcn.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h ${DIR}/photon-wave-function.o
	gcc ${FUNCDIR}/fcn.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${LIBS} -o ${DIR}/fcn ${INC} ${OPT} 
 
diffraction: ${DEPS} ${FUNCDIR}/diffraction.cc ${DIR}/dipole-cross-section.o ${FUNCDIR}/clenshaw.h ${FUNCDIR}/clenshaw-curtis.hh 
	g++ ${FUNCDIR}/diffraction.cc ${DIR}/dipole-cross-section.o ${LIBS} -o ${DIR}/diffraction ${INC}  ${OPT}

#ww-gluon-x: ${DEPS} ${UTILDIR}/ww-gluon-x.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/ww-gluon-x.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/ww-gluon-x   ${INC} ${OPT}

#ww-critical: ${DEPS} ${UTILDIR}/ww-critical.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/ww-critical.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/ww-critical   ${INC} ${OPT}

#geometric-scaling: ${DEPS} ${UTILDIR}/geometric-scaling.c  ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/geometric-scaling.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/geometric   ${INC} ${OPT}
#geometric: ${DEPS} ${UTILDIR}/geometric-scaling.c  ${UTILDIR}/tmd-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/geometric-scaling.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/geometric   ${INC} ${OPT}

#heavy-F2: ${DEPS} ${UTILDIR}/Plot-heavy-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o  
#	gcc ${UTILDIR}/Plot-heavy-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o ${LIBS} -o ${DIR}/heavy-F2 ${INC} ${OPT} 

export-data: ${DEPS} ${UTILDIR}/export-data.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${FUNCDIR}/clenshaw.h
	gcc ${UTILDIR}/export-data.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/export-data ${INC} ${OPT} 

#${DIR}/kt-formula.o: ${FUNCDIR}/kt-formula.c
#	gcc ${FUNCDIR}/kt-formula.c -o ${DIR}/kt-formula.o -c  ${INC} ${OPT} 


#${DIR}/clenshaw.o: ${FUNCDIR}/clenshaw.c
#	gcc -c ${FUNCDIR}/clenshaw.c -o ${DIR}/clenshaw.o  ${INC} ${OPT}


${DIR}/read-and-fit.o: ${FUNCDIR}/read-and-fit.c ${FUNCDIR}/read-and-fit-cheb.h ${FUNCDIR}/read-and-fit.h  ${UTILDIR}/f2.h ${FUNCDIR}/Kahn.h 
	gcc -c ${FUNCDIR}/read-and-fit.c -o ${DIR}/read-and-fit.o  ${INC} ${OPT}  

       
${DIR}/dipole-cross-section.o: ${FUNCDIR}/dipole-cross-section.c ${GLUDEP} ${SUDDEP}  
	gcc -c ${FUNCDIR}/dipole-cross-section.c -o ${DIR}/dipole-cross-section.o  ${INC} ${OPT}


${DIR}/photon-wave-function.o: ${FUNCDIR}/photon-wave-function.c 
	gcc -c ${FUNCDIR}/photon-wave-function.c  -o ${DIR}/photon-wave-function.o ${INC}  ${OPT}

		
clean: 
	rm -f ${DIR}/*.o
	rm -f ${DIR}/main
	 		
