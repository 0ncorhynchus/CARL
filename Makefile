#VPATH = src include
CPPFLAGS = -std=c++0x -g
LINK.o = g++

.PHONY: all
all: filter read_test

.PHONY: clean
clean:
	rm -f *.o filter read_test fasta_test

filter: read.o fasta.o filter.o filter_main.o
	g++ read.o fasta.o filter.o filter_main.o -lboost_thread-mt -o filter

read_test: read.o
fasta_test: read.o fasta.o

fasta.o: fasta.hpp
filter.o: filter.hpp
read.o: read.hpp
read_test.o:
fasta_test.o:
filter_main.o:
