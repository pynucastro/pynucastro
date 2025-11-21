CPPSOURCES := $(wildcard *.cpp)
FSOURCES := $(wildcard *.f90)

CPPOBJECTS := $(CPPSOURCES:.cpp=.o)
FOBJECTS := $(FSOURCES:.f90=.o)

HEADERS := $(wildcard *.H)

%.o: %.cpp $(HEADERS)
	g++ -I. -std=c++20 -c $<

%.o: %.f90
	gfortran -c $<

test: $(CPPOBJECTS) $(FOBJECTS)
	gfortran -o test $(FOBJECTS) $(CPPOBJECTS) -lstdc++


print-%:
	@echo "$* is $($*)"
