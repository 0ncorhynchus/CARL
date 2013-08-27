// fasta.hpp
// written by S.Kato

#ifndef __FASTA_hpp
#define __FASTA_hpp

#include <fstream>
#include <string>

typedef std::string Read;

class FastaItem {
private:
	std::string info;
	Read read;
public:
	FastaItem(std::string info, Read read);
	std::string getInfo();
	Read getRead();
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
