#VPATH = src include
CPPFLAGS = -std=c++0x -g
LINK.o = g++

.PHONY: all
all: filter

.PHONY: clean
clean:
	rm -f *.o filter read_test

filter: read.o fasta.o filter.o
	g++ read.o fasta.o filter.o -lboost_thread-mt -o filter

read_test: read.o read_test.o
	g++ read_test.o read.o -o read_test

fasta.o: fasta.hpp
filter.o: filter.hpp
read.o: read.hpp
read_test.o:

