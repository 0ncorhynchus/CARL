#include "read.hpp"
#include <random>
#include <iostream>

int main() {
	const char bases[5] = {'a', 'c', 'g', 't', 'n'};

	int length(13);
	char c_seq[length];
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> rnd(0, 4);
	for (int i(0); i < length; i++) {
		c_seq[i] = bases[rnd(mt)];
	}
	std::string seq(c_seq);
	std::cout << "====================" << std::endl;
	Read read(seq);
	std::cout << "ORIGIN: " << seq << std::endl;
	std::cout << "RESULT: " << read.tostring() <<  std::endl;
	bool isDefinite(read.isDefinite());
	std::cout << "IS DEFINITE: " << isDefinite << std::endl;

	std::cout << "====================" << std::endl;
	std::cout << "SUB: from 2 to 7" << std::endl;
	Read sub(read.sub(2, 6));
	std::cout << "ORIGIN: " << seq.substr(2, 6) << std::endl;
	std::cout << "RESULT: " << sub.tostring() << std::endl;

	std::cout << "====================" << std::endl;
	std::cout << "COMPLEMENT" << std::endl;
	Read comp(read.complement());
	std::cout << "ORIGIN: " << seq << std::endl;
	std::cout << "RESULT: " << comp.tostring() << std::endl;

	std::cout << "====================" << std::endl;
	std::string fseq("TCAGGGGGGTTTTAATTTACTTTCGTACACAGCGTAAATCTTACTAAATGTCTTACTATAACGCATACGATATCTTAACAACATCTAACTTCTAAAACATAGCACATTAAGCTCGAAAAACCAGCAAGCAAGCATACGAAGAAGTAAGAAATAATAACTCAATGTCGCTTCATTTTCTAGTTTAAAC");
	Read fasta(fseq);
	std::cout << fseq << std::endl;
	std::cout << fasta.tostring() << std::endl;

	std::cout << "====================" << std::endl;
	std::cout << "COPY" << std::endl;
	Read copy(fasta);
	std::cout << fasta.tostring() << std::endl;
	std::cout << copy.tostring() << std::endl;
}

