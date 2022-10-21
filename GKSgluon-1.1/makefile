# Makefile for GKSgluon-example

INCLUDE = $(shell gsl-config --cflags)
LIBRARIES+= $(shell gsl-config --libs-without-cblas)

CXX=g++
CXXFLAGS = -O -Wall -Wno-format -Wno-unused-variable -g -MMD -MP $(INCLUDE) 

# GSL
LIBRARIES+=-lgsl -lgslcblas

COMMONSRC = UnintegratedGluonGKS.cc InterpolationGKS.cc
COMMONOBJ = $(patsubst %.cc,%.o,$(COMMONSRC))
PROGSRC = GKSgluon-example.cc
PROGOBJ = $(patsubst %.cc,%.o,$(PROGSRC))
PROGEXE = $(patsubst %.cc,%,$(PROGSRC))
DEPS = $(patsubst %.cc,%.d,$(PROGSRC))

all: $(PROGEXE)

GKSgluon-example: $(COMMONOBJ) GKSgluon-example.o
	$(CXX) $(COMMONOBJ) -o $@ $@.o  $(LIBRARIES)

clean:
	rm -f $(COMMONOBJ) $(PROGOBJ) $(PROGEXE) *.d

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# takes care of dependencies
-include ${DEPS}
