SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
HEADERS := $(wildcard *.H)

%.o: %.cpp
	g++ -I. -c $<

main: $(OBJECTS) $(HEADERS)
	g++ -I. -o $@ $(OBJECTS)
