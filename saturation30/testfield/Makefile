	
all: integral simps.o


integral: simps.o integral.cpp simpson1dv1.hh
	g++ integral.cpp simps.o -o integral
	
#integral.o: simps.o 
#	g++ -c integral.cpp 
	

simps.o: ../simps.f
	gfortran -c ../simps.f

clean:	
	rm *.o
