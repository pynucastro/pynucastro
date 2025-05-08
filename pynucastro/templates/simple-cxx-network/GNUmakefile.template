SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
HEADERS := $(wildcard *.H)

DEFINES :=

ifeq ($(DISABLE_SCREENING),TRUE)
   DEFINES += -DDISABLE_SCREENING
endif

%.o: %.cpp
	g++ -I. $(DEFINES) -c $<

main: $(OBJECTS) $(HEADERS)
	g++ -I. -o $@ $(OBJECTS)
