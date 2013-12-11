// read.hpp
// written by S.Kato

#ifndef __READ_hpp
#define __READ_hpp

#include <string>
#include <stdexcept>

class Read {
private:
	static const char bases[4];

	unsigned char *_read;
	unsigned char *_flgs;
	const unsigned int _size;

	unsigned int r_size() const;
	unsigned int f_size() const;

public:
	Read(const std::string sequence);
	Read(const Read& read);
	~Read();
	unsigned int size() const;
	unsigned int getBaseAt(const unsigned int index) const throw(std::out_of_range);
	bool isDefinite() const;
	Read sub(const unsigned int start, const unsigned int length) const throw(std::out_of_range);
	Read complement() const;
	std::string tostring();

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
	Read(const unsigned int size);
	void setBaseAt(const unsigned int index, const unsigned int value) throw(std::out_of_range);
};

std::size_t hash_value(const Read& read);

#endif
