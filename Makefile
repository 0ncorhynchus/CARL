#VPATH = src include
CPPFLAGS = -std=c++0x
LINK.o = g++

.PHONY: all
all: filter

.PHONY: clean
clean:
	rm -f *.o filter

filter: fasta.o filter.o
	g++ fasta.o filter.o -lboost_thread-mt -o filter
fasta.o:
filter.o:
