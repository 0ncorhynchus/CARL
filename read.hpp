// read.hpp
// written by S.Kato

#ifndef __READ_hpp
#define __READ_hpp

#include <string>
#include <stdexcept>

class Read {
private:
	unsigned char *_read;
	unsigned char *_flgs;
	const unsigned int _size;
	static char bases[4];

public:
	Read(const std::string sequence);
	Read(const Read& read);
	~Read();
	unsigned int size() const;
	int getBaseAt(const unsigned int index) const throw(std::out_of_range);
	Read sub(const unsigned int start, const unsigned int length) const throw(std::out_of_range);
	Read complement() const;
	std::string tostring();

	bool operator==(const Read& read) {
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
	Read(const unsigned int size);
	void setBaseAt(const unsigned int index, const unsigned int value) throw(std::out_of_range);
};

char Read::bases[4] = {'a', 'c', 'g', 't'};

std::size_t hash_value(const Read& read) {
	std::size_t h = 1;
	for (int i(0); i < read.size(); i++) {
		h = h << 2;
		int base = read.getBaseAt(i);
		if (base & 4 == 1) {
			return 0;
		}
		h += (std::size_t)(base & 3);
	}
	return h;
}

#endif
