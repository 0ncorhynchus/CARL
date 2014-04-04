// fasta.cpp
// written by S.Kato

#include "fasta.hpp"
#include <regex.h>
#include <iostream>

namespace carl {

/*
 * Fasta::Item
 */
Fasta::Item::Item(std::string info, std::string sequence) : read(sequence) {
    this->info = info;
}

Fasta::Item::Item(const Fasta::Item& item) : read(item.read) {
    this->info = item.info;
}

Fasta::Item::Item() : read() {
    info = "";
}

Fasta::Item::~Item() {
}

std::string Fasta::Item::getInfo() const {
    return info;
}

const Read& Fasta::Item::getRead() const {
    return read;
}

/*
 * Fasta
 */
Fasta::Fasta(const std::string& filename) : _filename(filename){
    _ifs.open(_filename);
    getItemStrings();
}

Fasta::Fasta(const Fasta& fasta) : _filename(fasta._filename){
    _ifs.open(_filename);
    getItemStrings();
}

Fasta::~Fasta() {
    _ifs.close();
}

Fasta::Item Fasta::getItem() {
    std::pair<std::string, std::string> strings(this->getItemStrings());
    Fasta::Item item(strings.first, strings.second);
    return item;
}

std::pair<std::string, std::string> Fasta::getItemStrings() {
    std::string info, sequence;
    std::getline(_ifs, info);
    std::getline(_ifs, sequence);

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
                std::cerr << "Error in Fasta::getItemStrings()" << std::endl;
                std::cerr << "\"info\": " << info << std::endl;
                std::cerr << "\"so\": " << so << std::endl;
                std::cerr << "\"eo\": " << eo << std::endl;
                std::cerr << str << std::endl << std::endl;
            }
        }
    }
    regfree(&regexBuffer);

    std::pair<std::string, std::string> retval(_tmp);
    _tmp = std::pair<std::string, std::string>(info, sequence);
    return retval;
}

bool Fasta::eof() const {
    return _ifs.eof();
}

} // carl
