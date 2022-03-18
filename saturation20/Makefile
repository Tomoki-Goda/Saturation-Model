#------------------------------------------------------------------------------#
#                               Makefile: GBW Model                            #
#------------------------------------------------------------------------------#

CC = gcc $(INCLUDE)
FC = gfortran
LD = gfortran

#CUBADIR=$(HOME)/heptools/Cuba-4.2
CUBADIR=$(HOME)/Tools/Cuba-4.2
INCLUDE += -I$(CUBADIR)/include


DIRNAME = fit_result
# CERNLIB comes from global Ubuntu installation
LIBS   = -lpacklib -lmathlib -lkernlib #-lgsl -lgslcblas
#LIBS   = -lpacklib -lmathlib -lkernlib  -L$(CUBADIR)/lib -lcuba

FFLAGS = -c 
CFLAGS = -c -Wall -Df2cFortran -ggdb -Wno-unused-variable
LFLAGS = -Wl,-static -Wl,-dy -lnsl -lcrypt -ldl -lgfortran $(LIBS)

#          radmul0.o\

OBJECTS = main.o\
          simps.o\
          simps2d.o\
          bessels.o

HEADERS = control.h\
          main.h\
	  functions.h\
	  complex.h\
	  chebyshev.h\
	  chebyshev3.h\
	  gluons.h\
	  test.h\
          graph.h

TARGET  = main

all:    testcheb main graph
#all:    $(TARGET) test
#	cp $(TARGET) $(HOME)/bin
#	mkdir $(DIRNAME)
#	./main 
#	gnuplot dipol_cs.plt
#	gnuplot dipol_cs2.plt
#	gnuplot dipol_scaling.plt
#	gnuplot f_2_charm.plt
##	gnuplot diff_graph.plt
#	gnuplot diff_zeus94.plt
#	gnuplot diff_xp99.plt
#	gnuplot alpha_gluon.plt
#	gnuplot gluon.plt
#	gnuplot f2chf2.plt
#	gnuplot flq2.plt
#	 ./cpfigs
#	fixbb diff_zeus94.eps
#	fixbb diff_zeus99.eps
#	fixbb diff_xp99.eps
#	cp *.eps $(DIRNAME)
#	cp *.res $(DIRNAME)
#	cp *.dat $(DIRNAME)
#	cp *.plt $(DIRNAME)
#	rm *.tmp 


main.o: main.c $(HEADERS)
	$(CC) $(CFLAGS) main.c 

main:   main.o simps.o simps2d.o bessels.o
	$(LD) -o $(TARGET) $(OBJECTS) $(LFLAGS)

test.o: test.c $(HEADERS)
	$(CC) $(CFLAGS) test.c 

test:   test.o simps.o simps2d.o bessels.o
	$(LD) -o test test.o simps.o simps2d.o bessels.o $(LFLAGS)

testcheb.o: testcheb.c $(HEADERS)
	$(CC) $(CFLAGS) testcheb.c 

testcheb: testcheb.o simps.o simps2d.o bessels.o
	$(LD) -o testcheb testcheb.o simps.o simps2d.o bessels.o $(LFLAGS)

graph.o: graph.c $(HEADERS)
	$(CC) $(CFLAGS) graph.c 

graph: graph.o simps.o simps2d.o bessels.o
	$(LD) -o graph graph.o simps.o simps2d.o bessels.o $(LFLAGS)

simps.o: simps.f 
	$(FC) -c simps.f 

simps2d.o: simps2d.c
	$(CC) -c simps2d.c 
bessels.o: bessels.c
	$(CC) -c bessels.c 
clean:
	rm -f $(TARGET) main.o *.o test


################################################################################
#
# UWAGI
#
# 1 . Zmiana miejscem dolaczenia radmul0.o z bibliotekami cernlib podowuje blad,
#     kompilator mowi o podwojnej definicji  funkcji radmul.
# 2.  Mozna byloby pomyslec o przekompiniowaniu calego cernlibu tak aby 
#     radmul ktory tam powstanie byl poprawny.
# 3.  Inna mozliwoscia jest dolaczenie tylko plikow obiektowych z funkcjami
#     ktorych uzywamy w programie.
#
################################################################################




