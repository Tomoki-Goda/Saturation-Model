# This makefile can be used to build all or some of the sample
# programs. To build all of them, use the command
# 'make all'. To build one, say addtwo, use the command
# 'make addtwo'.

WSTPLINKDIR = /usr/local/Wolfram/Mathematica/12.0/SystemFiles/Links/WSTP/DeveloperKit
SYS = Linux-x86-64# Set this value with the result of evaluating $SystemID
CADDSDIR = ${WSTPLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR} 
FUNCDIR = ../Functions
LIBDIR = ${CADDSDIR}

EXTRALIBS = -lm -lpthread -lrt -lstdc++ -ldl -lgsl -lgslcblas -luuid # Set these with appropriate libs for your system.
WSTPLIB = WSTP64i4 # Set this to WSTP64i4 if using a 64-bit system

WSPREP = ${CADDSDIR}/wsprep
CXX=g++


all : gluons-math levin

gluons-math : gluons-math.o gluons.o clenshaw.o
	${CXX} gluons.o -I${INCDIR} gluons-math.o clenshaw.o -L${LIBDIR} -l${WSTPLIB} -L/usr/lib/x86_64-linux-gnu ${EXTRALIBS} -o $@

gluons.o: ${FUNCDIR}/gluons.cc ${FUNCDIR}/gluons.hh
	${CXX} -c ${FUNCDIR}/gluons.cc -o $@

clenshaw.o: ${FUNCDIR}/clenshaw2.cc ${FUNCDIR}/clenshaw2.hh
	${CXX} -c ${FUNCDIR}/clenshaw2.cc -o $@

levin : levin.o 
	${CXX} -I${INCDIR} levin.o -L${LIBDIR} -l${WSTPLIB} -L/usr/lib/x86_64-linux-gnu ${EXTRALIBS} -o $@

.c.o :
	${CXX} -c -I${INCDIR} -I${FUNCDIR}  $<

gluons-math.c : gluons-math.tm
	${WSPREP} $? -o $@

levin.c : levin.tm
	${WSPREP} $? -o $@

clean: 
	rm -f *.c *.o
