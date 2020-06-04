CXX=g++
# For Linux
CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -march=native -fopenmp
  
.DEFAULT_GOAL := all

all: Lock 
 
TAS_lock: Lock.cpp Makefile
	$(CXX) Lock.cpp -o Lock -v  $(CXXFLAGS)  

.PHONY: clean
clean:
	rm Lock
