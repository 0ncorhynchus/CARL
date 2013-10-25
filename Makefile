#VPATH = src include
CPPFLAGS = -std=c++0x -g
LINK.o = g++

.PHONY: all
all: filter

.PHONY: clean
clean:
	rm -f *.o filter read

filter: fasta.o filter.o
	g++ fasta.o filter.o -lboost_thread-mt -o filter

read:
fasta.o:
filter.o:
read.o:

