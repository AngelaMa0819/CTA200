XX = g++
CXXFLAGS = -Wall -O2
LDFLAGS = -lm

all: average

average: average.o
        $(CXX) -o $@ $< $(LDFLAGS)

%.o: %.cc
        $(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
        $(RM) average *.o

