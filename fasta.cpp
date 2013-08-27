// fasta.cpp
// written by S.Kato

#include "fasta.hpp"
#include <regex.h>

Fasta::Fasta(std::string& filename) : ifs(filename.c_str()) {
}

Fasta::~Fasta() {
	ifs.close();
}

FastaItem Fasta::getItem() {
	std::string info, sequence;
	std::getline(ifs, info);
	std::getline(ifs, sequence);

	const char regex[] = "^>\\b(.*)";
	regex_t regexBuffer;
	if (regcomp(&regexBuffer, regex, REG_EXTENDED | REG_NEWLINE) == 0) {
		regmatch_t patternMatch[2];
		int size = sizeof(patternMatch) / sizeof(regmatch_t);
		if (regexec(&regexBuffer, info.c_str(), size, patternMatch, 0) == 0) {
			int so = patternMatch[1].rm_so;
			int eo = patternMatch[1].rm_eo;
			if (so != -1 && eo != -1) {
				info = info.substr(so, eo-so);
			}
		}
	}

	FastaItem item(info, sequence);
	regfree(&regexBuffer);
	return item;
}

bool Fasta::eof() {
	return ifs.eof();
}

FastaItem::FastaItem(std::string info, Read read) {
	this->info = info;
	this->read = read;
}

std::string FastaItem::getInfo() {
	return info;
}

Read FastaItem::getRead() {
	return read;
}
