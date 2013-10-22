// fasta.cpp
// written by S.Kato

#include "fasta.hpp"
#include <regex.h>
#include <iostream>

Fasta::Fasta(std::string& filename) : ifs(filename.c_str()) {
}

Fasta::~Fasta() {
	ifs.close();
}


FastaItem Fasta::getItem() {
	std::lock_guard<std::mutex> lock(m_mutex);
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
			try {
				info = info.substr(so, eo-so);
			} catch(const char* str) {
				std::cerr << "Error in Fasta::getItem()" << std::endl;
				std::cerr << str << std::endl;
			}
		}
	}

	FastaItem item(info, sequence);
	regfree(&regexBuffer);
	return item;
}

bool Fasta::eof() {
	std::lock_guard<std::mutex> lock(m_mutex);
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
