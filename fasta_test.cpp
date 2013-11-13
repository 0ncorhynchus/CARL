#include "fasta.hpp"
#include <string>
#include <iostream>

int main() {
	std::string filename("test.fasta");
	Fasta fasta(filename);
	while (!fasta.eof()) {
		FastaItem item(fasta.getItem());
		std::cout << "> " << item.getInfo() << std::endl;
		Read read(item.getRead());
		std::cout << read.tostring() << std::endl;
	}
}
