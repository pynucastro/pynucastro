SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
HEADERS := $(wildcard *.H)

DEFINES :=

ifeq ($(DISABLE_SCREENING),TRUE)
   DEFINES += -DDISABLE_SCREENING
endif

%.o: %.cpp $(HEADERS)
	g++ -std=c++17 -I. $(DEFINES) -c $<

main: $(OBJECTS) $(HEADERS)
	g++ -std=c++17 -I. -o $@ $(OBJECTS)
