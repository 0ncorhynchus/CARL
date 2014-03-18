// read.hpp
// written by S.Kato

#ifndef __READ_hpp
#define __READ_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>

class Read {
private:
	static const char bases[4];

	std::vector<unsigned char> _read;
	std::vector<unsigned char> _flgs;

	unsigned int _size;

public:
	Read(const std::string sequence);
	Read(const Read& read);
	Read(const unsigned int size);
	Read();
	unsigned int size() const;
	unsigned char getBaseAt(const unsigned int index) const throw(std::out_of_range);
	bool isDefinite() const;
	Read sub(const unsigned int start, const unsigned int length) const
		throw(std::out_of_range);
	Read complement() const;
	Read reverse() const;
	std::string tostring() const;

	bool operator==(const Read& read) const {
		if (this->size() != read.size())
			return false;
		bool flg(true);
		for (int i(0); i < this->size(); i++) {
			if (this->getBaseAt(i) != read.getBaseAt(i)) {
				flg = false;
				break;
			}
		}
		return flg;
	}

private:
	std::pair<unsigned int, unsigned int> _indexes(unsigned int index) const;
	void setBaseAt(const unsigned int index, const unsigned char value)
		throw(std::out_of_range);
};

std::size_t hash_value(const Read& read);

#endif
