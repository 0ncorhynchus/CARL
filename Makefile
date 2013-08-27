#VPATH = src include
#CPPFLAGS = -I include
LINK.o = g++

.PHONY: all
all: filter

.PHONY: clean
clean:
	rm -f *.o filter

filter: fasta.o filter.o
fasta.o:
filter.o:
