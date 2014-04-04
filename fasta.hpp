// fasta.hpp
// written by S.Kato

#ifndef __FASTA_hpp
#define __FASTA_hpp

#include <fstream>
#include <string>
#include "read.hpp"

namespace carl {

class Fasta {
public:
    class Item {
        private:
            std::string info;
            Read read;
        public:
            Item(std::string info, std::string sequence);
            Item(const Item& item);
            Item();
            ~Item();
            std::string getInfo() const;
            const Read& getRead() const;
    };

private:
    const std::string _filename;
    std::ifstream _ifs;
    std::pair<std::string, std::string> _tmp;
public:
    Fasta(const std::string& filename);
    Fasta(const Fasta& fasta);
    ~Fasta();
    Item getItem();
    std::pair<std::string, std::string> getItemStrings();
    bool eof() const;
};

} // carl

#endif
