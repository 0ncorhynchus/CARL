#include "read.hpp"

namespace carl {

const char Read::bases[4] = {'a', 'c', 'g', 't'};

Read::Read(const std::string sequence) : _size(sequence.size()) {
    unsigned char read(0), flg(0);
    for (size_type i(0); i < size(); i++)
    {
        unsigned int bp(4);
        switch(sequence.at(i)) {
            case 'a':
            case 'A':
                bp = 0;
                break;

            case 'c':
            case 'C':
                bp = 1;
                break;

            case 'g':
            case 'G':
                bp = 2;
                break;

            case 't':
            case 'T':
                bp = 3;
                break;
        }
        read += (bp & 3) << ((i & 3) * 2);
        if ((i & 3) == 3) {
            _read.push_back(read);
            read = 0;
        }
        flg += ((bp >> 2) & 1) << (i & 7);
        if ((i & 7) == 7) {
            _flgs.push_back(flg);
            flg = 0;
        }
    }
    if ((size() & 3) != 0) {
        _read.push_back(read);
    }
    if ((size() & 7) != 0) {
        _flgs.push_back(flg);
    }
}

Read::Read(const Read& read) : _size(read.size()) {
    _read.resize(read._read.size());
    _flgs.resize(read._flgs.size());
    std::copy(read._read.begin(), read._read.end(), _read.begin());
    std::copy(read._flgs.begin(), read._flgs.end(), _flgs.begin());
}

Read::Read(const size_type size) : _size(size) {
    const std::pair<unsigned int, unsigned int> sizes(_indexes(size));
    _read.resize(sizes.first+1, 0);
    _flgs.resize(sizes.second+1, 0);
}

Read::Read() : _size(0) {
}

std::pair<Read::size_type, Read::size_type> Read::_indexes(size_type index) const {
    std::pair<size_type, size_type> retval(0,0);
    retval.first = index >> 2;
    retval.second = index >> 3;
    return retval;
}

Read::size_type Read::size() const {
    return this->_size;
}

unsigned char Read::getBaseAt(const size_type index) const
throw(std::out_of_range) {
    if (index >= this->size())
        throw std::out_of_range("out of range in getBaseAt()");

    unsigned char score(0);

    const int f_index(index >> 3);
    const int f_order(index & 7);
    score += ((_flgs.at(f_index) >> f_order) & 1) << 2;

    const int r_index(index >> 2);
    const int r_order((index & 3) << 1);
    score += (_read.at(r_index) >> r_order) & 3;
    return score;
}

void Read::setBaseAt(const size_type index, const unsigned char value)
throw(std::out_of_range) {
    if (index >= this->size())
        throw std::out_of_range("out of range in setBaseAt()");

    const int r_index(index >> 2);
    const int r_order((index & 3) << 1);
    unsigned char read(_read.at(r_index));
    read &= 255 - (3 << r_order);
    read += (value & 3) << r_order;
    _read[r_index] = read;

    const int f_index(index >> 3);
    const int f_order(index & 7);
    unsigned char flg(_flgs.at(f_index));
    flg &= 255 - (1 << f_order);
    flg += ((value >> 2) & 1) << f_order;
    _flgs[f_index] = flg;
}

bool Read::isDefinite() const {
    for (std::vector<unsigned char>::const_iterator itr(_flgs.begin());
            itr != _flgs.end(); itr++) {
        if ((unsigned int)(*itr) != 0) {
            return false;
        }
    }
    return true;
}

Read Read::sub(const size_type start, const size_type length) const throw(std::out_of_range){
    if (length <= 0 || start + length > this->size())
        throw std::out_of_range("out of range in sub()");

    Read retval(length);
    for (size_type i(0); i < length; i++) {
        int base(this->getBaseAt(start + i));
        retval.setBaseAt(i, base);
    }
    return retval;
}

Read Read::complement() const {
    if (this->size() == 0)
        return Read();
    Read retval(*this);
    for (std::vector<unsigned char>::iterator itr(retval._read.begin());
            itr != retval._read.end(); itr++) {
        (*itr) = ~(*itr);
    }
    return retval;
}

Read Read::reverse() const {
    if (size() == 0)
        return Read();
    Read retval(*this);
    for (size_type i(0); i < size(); i++) {
        unsigned char base(this->getBaseAt(i));
        retval.setBaseAt(size()-i-1, base);
    }
    return retval;
}

std::string Read::tostring() const {
    std::string str("");

    for (size_type i(0); i < size(); i++) {
        int bp(this->getBaseAt(i));
        char ch('n');
        if (bp >= 0 && bp < 4)
            ch = Read::bases[bp];
        str += ch;
    }

    return str;
}

} // carl

std::size_t hash_value(const carl::Read& read) {
    std::size_t h(1);
    for (carl::Read::size_type i(0); i < read.size(); i++) {
        h = h << 2;
        int base = read.getBaseAt(i);
        if ((base & 4) == 1) {
            return 0;
        }
        h += (std::size_t)(base & 3);
    }
    return h;
}

