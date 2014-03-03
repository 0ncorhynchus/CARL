// fasta.cpp
// written by S.Kato

#include "fasta.hpp"
#include <regex.h>
#include <iostream>

/*
 * FastaItem
 */
FastaItem::FastaItem(std::string info, std::string sequence) : read(sequence) {
	this->info = info;
}

FastaItem::FastaItem(const FastaItem& item) : read(item.read) {
	this->info = item.info;
}

FastaItem::~FastaItem() {
}

std::string FastaItem::getInfo() const {
	return info;
}

const Read& FastaItem::getRead() const {
	return read;
}

/*
 * Fasta
 */
Fasta::Fasta(const std::string& filename) : filename(filename){
	this->ifs.open(this->filename);
}

Fasta::Fasta(const Fasta& fasta) : filename(fasta.filename){
	ifs.open(this->filename);
}

Fasta::~Fasta() {
	ifs.close();
}


FastaItem Fasta::getItem() {
	std::pair<std::string, std::string> strings(this->getItemStrings());
	FastaItem item(strings.first, strings.second);
	return item;
}

std::pair<std::string, std::string> Fasta::getItemStrings() {
	std::string info, sequence;
	std::getline(ifs, info);
	std::getline(ifs, sequence);

	const char regex[] = "^>(.*)";
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
				std::cerr << "\"info\": " << info << std::endl;
				std::cerr << "\"so\": " << so << std::endl;
				std::cerr << "\"eo\": " << eo << std::endl;
				std::cerr << str << std::endl << std::endl; 
			}
		}
	}

	std::pair<std::string, std::string> retval(info, sequence);
	regfree(&regexBuffer);
	return retval;
}

bool Fasta::eof() const {
	return ifs.eof();
}
