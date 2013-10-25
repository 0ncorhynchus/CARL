#include "read.hpp"

inline const unsigned int getReadIndex(unsigned int index) {
	return index >> 2;
}

inline const unsigned int getFlgsIndex(unsigned int index) {
	return index >> 3;
}

Read::Read(const std::string sequence) : _size(sequence.size()) {
	if (this->size() == 0) {
		this->_read = new unsigned char[0];
		this->_flgs = new unsigned char[0];
		return;
	}

	this->_read = new unsigned char[((this->size() - 1) >> 2) + 1];
	this->_flgs = new unsigned char[((this->size() - 1) >> 3) + 1];

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
	int r_length(sizeof(read._read));
	int f_length(sizeof(read._flgs));
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
		this->_read = new unsigned char[((this->size() - 1) >> 2) + 1];
		this->_flgs = new unsigned char[((this->size() - 1) >> 3) + 1];
	}
}

Read::~Read() {
	delete[] this->_read;
	delete[] this->_flgs;
}

unsigned int Read::size() const {
	return this->_size;
}

int Read::getBaseAt(const unsigned int index) const throw(std::out_of_range) {
	if (index >= this->size())
		throw std::out_of_range("out of range in getBaseAt()");

	const int f_index(index >> 3);
	const int f_order(index & 7);
	if ((this->_flgs[f_index] >> f_order) & 1 == 1)
		return -1;

	const int r_index(index >> 2);
	const int r_order((index & 3) * 2);
	return (this->_read[r_index] >> r_order) & 3;
}

void Read::setBaseAt(const unsigned int index, const unsigned int value) throw(std::out_of_range) {
	if (index >= this->size())
		throw std::out_of_range("out of range in setBaseAt()");

	const int r_index(index >> 2);
	const int r_order((index & 3) * 2);
	this->_read[r_index] &= 255 - (3 << r_order);
	this->_read[r_index] += (value & 3) << r_order;

	const int f_index(index >> 3);
	const int f_order(index & 7);
	this->_flgs[f_index] &= 255 - (1 << f_order);
	this->_flgs[f_index] += ((value >> 2) & 1) << f_order;
}

Read Read::sub(const unsigned int start, const unsigned int length) const throw(std::out_of_range){
	if (start < 0 || length <= 0 || start + length >= this->size())
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
	for (int i(0); i < sizeof(this->_read); i++) {
		retval._read[i] = ~this->_read[i];
	}
	for (int i(0); i < sizeof(this->_flgs); i++) {
		retval._flgs[i] = this->_flgs[i];
	}
	return retval;
}

std::string Read::tostring() {
	if (this->size() == 0)
		return std::string("");
	char* c_str = new char[this->size()];
	for (int i(0); i < this->size(); i++) {
		int bp(this->getBaseAt(i));
		char ch('n');
		if (bp >= 0 && bp < 4)
			ch = Read::bases[bp];
		c_str[i] = ch;
	}
	std::string str(c_str);
	delete[] c_str;
	return str;
}

#include <random>
#include <iostream>
int main() {
	char bases[5] = {'a', 'c', 'g', 't', 'n'};

	int length(13);
	char c_seq[length];
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_int_distribution<int> rnd(0, 4);
	for (int i(0); i < length; i++) {
		c_seq[i] = bases[rnd(mt)];
	}
	std::string seq(c_seq);
	std::cout << "====================" << std::endl;
	Read read(seq);
	std::cout << "ORIGIN: " << seq << std::endl;
	std::cout << "RESULT: " << read.tostring() <<  std::endl;

	std::cout << "====================" << std::endl;
	std::cout << "SUB: from 2 to 7" << std::endl;
	Read sub(read.sub(2, 6));
	std::cout << "ORIGIN: " << seq.substr(2, 6) << std::endl;
	std::cout << "RESULT: " << sub.tostring() << std::endl;

	std::cout << "====================" << std::endl;
	std::cout << "COMPLEMENT" << std::endl;
	Read comp(read.complement());
	std::cout << "ORIGIN: " << seq << std::endl;
	std::cout << "RESULT: " << comp.tostring() << std::endl;

}

