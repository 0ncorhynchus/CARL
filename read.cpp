#include "read.hpp"

const char Read::bases[4] = {'a', 'c', 'g', 't'};

std::size_t hash_value(const Read& read) {
	std::size_t h(1);
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

Read::Read(const std::string sequence) : _size(sequence.size()) {
	if (this->size() == 0) {
		this->_read = new unsigned char[0];
		this->_flgs = new unsigned char[0];
		return;
	}

	this->_read = new unsigned char[this->r_size()];
	this->_flgs = new unsigned char[this->f_size()];
	for (int i(0); i < this->r_size(); i++) {
		this->_read[i] = 0;
	}
	for (int i(0); i < this->f_size(); i++) {
		this->_flgs[i] = 0;
	}

	for (int i(0); i < this->size(); i++)
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

		this->setBaseAt(i, bp);
	}
}

Read::Read(const Read& read) : _size(read.size()) {
	int r_length(read.r_size());
	int f_length(read.f_size());

	this->_read = new unsigned char[r_length];
	this->_flgs = new unsigned char[f_length];

	for (int i(0); i < r_length; i++) {
		this->_read[i] = read._read[i];
	}
	for (int i(0); i < f_length; i++) {
		this->_flgs[i] = read._flgs[i];
	}
}

Read::Read(const unsigned int size) : _size(size) {
	if (this->size() == 0) {
		this->_read = new unsigned char[0];
		this->_flgs = new unsigned char[0];
	} else {
		int r_length(this->r_size());
		int f_length(this->f_size());
		this->_read = new unsigned char[r_length];
		this->_flgs = new unsigned char[f_length];
		for (int i(0); i < r_length; i++) {
			this->_read[i] = 0;
		}
		for (int i(0); i < f_length; i++) {
			this->_flgs[i] = 0;
		}
	}
}

Read::~Read() {
	delete[] this->_read;
	delete[] this->_flgs;
}

unsigned int Read::r_size() const {
	if (this->size() == 0) {
		return 0;
	}
	int size(((this->size() - 1) >> 2) + 1);
	return size;
}

unsigned int Read::f_size() const {
	if (this->size() == 0) {
		return 0;
	}
	int size(((this->size() - 1) >> 3) + 1);
	return size;
}

unsigned int Read::size() const {
	return this->_size;
}

unsigned int Read::getBaseAt(const unsigned int index) const throw(std::out_of_range) {
	if (index >= this->size())
		throw std::out_of_range("out of range in getBaseAt()");
	int score(0);

	const int f_index(index >> 3);
	const int f_order(index & 7);
	score += ((this->_flgs[f_index] >> f_order) & 1) << 2;

	const int r_index(index >> 2);
	const int r_order((index & 3) << 1);
	score += (this->_read[r_index] >> r_order) & 3;
	return score;
}

void Read::setBaseAt(const unsigned int index, const unsigned int value) throw(std::out_of_range) {
	if (index >= this->size())
		throw std::out_of_range("out of range in setBaseAt()");

	const int r_index(index >> 2);
	const int r_order((index & 3) << 1);
	this->_read[r_index] &= 255 - (3 << r_order);
	this->_read[r_index] += (value & 3) << r_order;

	const int f_index(index >> 3);
	const int f_order(index & 7);
	this->_flgs[f_index] &= 255 - (1 << f_order);
	this->_flgs[f_index] += ((value >> 2) & 1) << f_order;
}

bool Read::isDefinite() const {
	bool flg(true);
	for (int i(0); i < this->f_size(); i++) {
		if ((unsigned int)(this->_flgs[i]) != 0) {
			flg = false;
			break;
		}
	}
	return flg;
}

Read Read::sub(const unsigned int start, const unsigned int length) const throw(std::out_of_range){
	if (start < 0 || length <= 0 || start + length > this->size())
		throw std::out_of_range("out of range in sub()");

	Read retval(length);
	for (int i(0); i < length; i++) {
		int base(this->getBaseAt(start + i));
		retval.setBaseAt(i, base);
	}
	return retval;
}

Read Read::complement() const {
	if (this->size() == 0)
		return Read(0);

	Read retval(this->size());
	for (int i(0); i < this->size(); i++) {
		unsigned int base(this->getBaseAt(i));
		retval.setBaseAt(this->size()-i-1, base);
	}

	for (int i(0); i < retval.r_size(); i++) {
		retval._read[i] = ~retval._read[i];
	}
	return retval;
}

std::string Read::tostring() {
	int size(this->size());
	std::string str("");

	for (int i(0); i < size; i++) {
		int bp(this->getBaseAt(i));
		char ch('n');
		if (bp >= 0 && bp < 4)
			ch = Read::bases[bp];
		str += ch;
	}

	return str;
}

