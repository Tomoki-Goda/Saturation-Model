#########################################################
# before use do export DIR=<directory >			#
# like 							#
# export DIR						#
# for i in ./GBW/* ;  DIR=${i}; make ;done	        #
#							#
#########################################################


########################################################################################################
DEPS= ${DIR}/control.h  ./Functions/control-default.h ./Functions/constants.h ./Functions/Parameters.h 

LIBS=-lm -lmathlib -lkernlib -lpacklib  -Wno-unused-result -lquadmath

UTILDIR=./Utilities
FUNCDIR=./Functions

INC=-I ${DIR} -I ${FUNCDIR} -I ${UTILDIR}

OPT=-O3

SUDDEP=./Functions/mu2.h ./Functions/sudakov.h
GLUDEP=./Functions/gluon-chebyshev.h ./Functions/gluons.h ./Functions/chebyshev-1.h ./Functions/complex.h

########################################################################################################

all: clean main #critical dipole F2 FL F2-slope psi Integrand tmd-gluon tmd-gluon-x tmd-critical fcn export-data geometric ww-gluon ww-gluon-x ww-critical heavy-F2 
grid: gluon-grid

fit: clean main

plot: critical dipole F2 FL F2-slope tmd-gluon tmd-gluon-x tmd-critical fcn export-data  ww-gluon 


main: ${DEPS} ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o ${DIR}/clenshaw.o 
	gcc ${FUNCDIR}/main.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o ${LIBS} -o ${DIR}/main ${INC} ${OPT} 


gluon-grid: ${DEPS} ${UTILDIR}/gluon-grid.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon.h
	gcc ${UTILDIR}/gluon-grid.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/gluon-grid   ${INC} ${OPT}


dipole-grid: ${DEPS} ${UTILDIR}/dipole-grid.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o
	gcc ${UTILDIR}/dipole-grid.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/dipole-grid  ${INC} ${OPT}

critical : ${DEPS} ${UTILDIR}/critical-line.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o
	gcc ${UTILDIR}/critical-line.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/critical   ${INC} ${OPT}

dipole: ${DEPS} ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o
	gcc ${UTILDIR}/Plot-DP.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/dipole  ${INC} ${OPT}

F2: ${DEPS} ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o ${DIR}/clenshaw.o 
	gcc ${UTILDIR}/Plot-F2.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o ${LIBS} -o ${DIR}/F2   ${INC} ${OPT}

FL: ${DEPS} ${UTILDIR}/Plot-FL.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o ${DIR}/clenshaw.o 
	gcc ${UTILDIR}/Plot-FL.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/FL   ${INC} ${OPT}

F2-slope: ${DEPS} ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o ${DIR}/clenshaw.o 
	gcc ${UTILDIR}/F2-slope.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/F2-slope   ${INC} ${OPT}

#psi: ${DEPS} ${UTILDIR}/Plot-Psi.c ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o 
#	gcc ${UTILDIR}/Plot-Psi.c  ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/psi   ${INC} ${OPT}

curvature : ${DEPS} ${UTILDIR}/curvature.c ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o
	gcc ${UTILDIR}/curvature.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/curvature  ${INC} ${OPT}

#Integrand: ${DEPS} ${UTILDIR}/Plot-Integrand.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o 
#	gcc ${UTILDIR}/Plot-Integrand.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${LIBS} -o ${DIR}/Integrand   ${INC} ${OPT}

tmd-gluon: ${DEPS} ${UTILDIR}/tmd-gluon.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon.h ${DIR}/clenshaw.o
	gcc ${UTILDIR}/tmd-gluon.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/tmd-gluon   ${INC} ${OPT}

tmd-gluon-x: ${DEPS} ${UTILDIR}/tmd-gluon-x.c ${DIR}/dipole-cross-section.o  ${UTILDIR}/tmd-gluon.h ${DIR}/clenshaw.o
	gcc ${UTILDIR}/tmd-gluon-x.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/tmd-gluon-x   ${INC} ${OPT}

tmd-critical: ${DEPS} ${UTILDIR}/tmd-critical.c  ${UTILDIR}/tmd-gluon.h ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o 
	gcc ${UTILDIR}/tmd-critical.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/tmd-critical   ${INC} ${OPT}

ww-gluon: ${DEPS} ${UTILDIR}/ww-gluon.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o 
	gcc ${UTILDIR}/ww-gluon.c ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o  ${LIBS} -o ${DIR}/ww-gluon   ${INC} ${OPT}

#ww-gluon-x: ${DEPS} ${UTILDIR}/ww-gluon-x.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/ww-gluon-x.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/ww-gluon-x   ${INC} ${OPT}

#ww-critical: ${DEPS} ${UTILDIR}/ww-critical.c ${UTILDIR}/ww-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/ww-critical.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/ww-critical   ${INC} ${OPT}

#geometric-scaling: ${DEPS} ${UTILDIR}/geometric-scaling.c  ${UTILDIR}/critical-line.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/geometric-scaling.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/geometric   ${INC} ${OPT}
#geometric: ${DEPS} ${UTILDIR}/geometric-scaling.c  ${UTILDIR}/tmd-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/geometric-scaling.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/geometric   ${INC} ${OPT}

#tmd-critical-2: ${DEPS} ${UTILDIR}/tmd-critical-2.c  ${UTILDIR}/tmd-gluon.h ${DIR}/dipole-cross-section.o 
#	gcc ${UTILDIR}/tmd-critical-2.c ${DIR}/dipole-cross-section.o  ${LIBS} -o ${DIR}/tmd-critical-2   ${INC} ${OPT}

fcn: ${DEPS} ${FUNCDIR}/fcn.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/clenshaw.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o ${DIR}/dcurtis.o 
	gcc ${FUNCDIR}/fcn.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/clenshaw.o ${DIR}/photon-wave-function.o   ${DIR}/dcurtis.o ${LIBS} -o ${DIR}/fcn ${INC} ${OPT} 

export-data: ${DEPS} ${UTILDIR}/export-data.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o ${DIR}/simpson-integral.o  ${DIR}/dcurtis.o
	gcc ${UTILDIR}/export-data.c ${DIR}/read-and-fit.o ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o ${DIR}/clenshaw.o  ${DIR}/dcurtis.o ${LIBS} -o ${DIR}/export-data ${INC} ${OPT} 
 
heavy-F2: ${DEPS} ${UTILDIR}/Plot-heavy-F2.c ${DIR}/dipole-cross-section.o ${DIR}/photon-wave-function.o ${DIR}/simpson-integral.o  ${DIR}/dcurtis.o ${DIR}/clenshaw.o  
	gcc ${UTILDIR}/Plot-heavy-F2.c ${DIR}/dipole-cross-section.o ${DIR}/simpson-integral.o ${DIR}/photon-wave-function.o  ${DIR}/dcurtis.o ${DIR}/clenshaw.o ${LIBS} -o ${DIR}/heavy-F2 ${INC} ${OPT} 


${DIR}/dgauss.o: ${FUNCDIR}/dgauss.c
	gcc -c ${FUNCDIR}/dgauss.c ${LIBS} ${INC} -o ${DIR}/dgauss.o ${OPT}

${DIR}/dcurtis.o: ${FUNCDIR}/dcurtis.c
	gcc -c ${FUNCDIR}/dcurtis.c ${LIBS} ${INC} -o ${DIR}/dcurtis.o ${OPT}

${DIR}/clenshaw.o: ${FUNCDIR}/clenshaw.c
	gcc -c ${FUNCDIR}/clenshaw.c ${LIBS} ${INC} -o ${DIR}/clenshaw.o ${OPT}

#${DIR}/clenshaw-curtis.o: ${FUNCDIR}/clenshaw-curtis.cc
#	g++ -c ${FUNCDIR}/clenshaw-curtis.cc ${LIBS} ${INC} -o ${DIR}/clenshaw-curtis.o 


diffraction: ${DEPS} ${FUNCDIR}/diffraction.cc ${DIR}/dipole-cross-section.o  ${DIR}/clenshaw.o  ${DIR}/dcurtis.o  ${DIR}/dgauss.o ${FUNCDIR}/clenshaw-curtis.hh 
	g++ ${FUNCDIR}/diffraction.cc ${DIR}/dipole-cross-section.o  ${DIR}/dgauss.o ${DIR}/clenshaw.o ${DIR}/dcurtis.o ${LIBS} -o ${DIR}/diffraction ${INC} -std=c++17 -O1


${DIR}/read-and-fit.o: ${DEPS} ${FUNCDIR}/read-and-fit.c ${UTILDIR}/f2.h
	echo ${DIR}
	gcc -c ${FUNCDIR}/read-and-fit.c   ${LIBS}  ${INC}  -o ${DIR}/read-and-fit.o ${OPT}

       
${DIR}/dipole-cross-section.o: ${DEPS} ${FUNCDIR}/dipole-cross-section.c ${GLUDEP} ${SUDDEP}  
	gcc -c ${FUNCDIR}/dipole-cross-section.c  ${LIBS}   ${INC}   -o ${DIR}/dipole-cross-section.o ${OPT}


${DIR}/photon-wave-function.o: ${DEPS}  ${FUNCDIR}/photon-wave-function.c 
	gcc -c ${FUNCDIR}/photon-wave-function.c  ${LIBS}   ${INC}  -o ${DIR}/photon-wave-function.o  ${OPT}


${DIR}/simpson-integral.o: ${FUNCDIR}/simpson-integral.c
	gcc -c ${FUNCDIR}/simpson-integral.c  ${LIBS}   ${INC} -o ${DIR}/simpson-integral.o  ${OPT}

		
clean: 
	rm -f ${DIR}/*.o
	rm -f ${DIR}/main
	 		
