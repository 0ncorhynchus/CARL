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
	const std::string filename;
	std::ifstream ifs;
public:
	Fasta(const std::string& filename);
	Fasta(const Fasta& fasta);
	~Fasta();
	FastaItem getItem();
	std::pair<std::string, std::string> getItemStrings();
	bool eof() const;
};

#endif
