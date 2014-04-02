// read.hpp
// written by S.Kato

#ifndef __READ_hpp
#define __READ_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>

class Read {
public:
    typedef unsigned int size_type;
private:
    static const char bases[4];

    std::vector<unsigned char> _read;
    std::vector<unsigned char> _flgs;

    size_type _size;

public:
    Read(const std::string sequence);
    Read(const Read& read);
    Read(const size_type size);
    Read();
    size_type size() const;
    unsigned char getBaseAt(const size_type index) const throw(std::out_of_range);
    bool isDefinite() const;
    Read sub(const size_type start, const size_type length) const
        throw(std::out_of_range);
    Read complement() const;
    Read reverse() const;
    std::string tostring() const;

    bool operator==(const Read& read) const {
        if (this->size() != read.size())
            return false;
        bool flg(true);
        for (size_type i(0); i < this->size(); i++) {
            if (this->getBaseAt(i) != read.getBaseAt(i)) {
                flg = false;
                break;
            }
        }
        return flg;
    }

private:
    std::pair<size_type, size_type> _indexes(size_type index) const;
    void setBaseAt(const size_type index, const unsigned char value)
        throw(std::out_of_range);
};

std::size_t hash_value(const Read& read);

namespace std {
    template <>
    struct hash<Read> {
        std::size_t operator()(const Read& read) const {
            return hash_value(read);
        }
    };
}

#endif
