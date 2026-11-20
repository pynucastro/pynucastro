SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
HEADERS := $(wildcard *.H)

DEFINES :=

USE_SCREENING ?= TRUE
ifeq ($(USE_SCREENING),TRUE)
   DEFINES += -DSCREENING
endif

%.o: %.cpp $(HEADERS)
	g++ -std=c++20 -I. $(DEFINES) -c $<

main: $(OBJECTS) $(HEADERS)
	g++ -std=c++20 -I. -o $@ $(OBJECTS)
