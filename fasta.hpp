// fasta.hpp
// written by S.Kato

#ifndef __FASTA_hpp
#define __FASTA_hpp

#include <fstream>
#include <string>
#include "read.hpp"

class FastaItem {
private:
	std::string info;
	Read read;
public:
	FastaItem(std::string info, std::string sequence);
	FastaItem(const FastaItem& item);
	~FastaItem();
	std::string getInfo() const;
	const Read& getRead() const;
};

class Fasta {
private:
	std::ifstream ifs;
public:
	Fasta(std::string& filename);
	~Fasta();
	FastaItem getItem();
	bool eof();
};

#endif
