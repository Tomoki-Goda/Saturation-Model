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


all : gluons-math 

gluons-math : gluons-math.o
	${CXX} -I${INCDIR} gluons-math.o -L${LIBDIR} -l${WSTPLIB} -L/usr/lib/x86_64-linux-gnu ${EXTRALIBS} -o $@

.c.o :
	${CXX} -c -I${INCDIR} -I${FUNCDIR}  $<

gluons-math.c : gluons-math.tm
	${WSPREP} $? -o $@
