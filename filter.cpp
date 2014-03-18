// filter.cpp
// written by S.Kato

#include <iostream>
#include <string>
#include "filter.hpp"

Filter::Filter(int lower_level, int lower_interval, int upper_level, double ratio) {
	this->_mer_length = 0;
	this->_lower_level = lower_level;
	this->_lower_interval = lower_interval;
	this->_upper_level = upper_level;
	this->_ratio = ratio;
	this->_debug = false;
}

Filter::Filter(const Filter& filter) {
	this->_mer_length = filter._mer_length;
	this->_lower_level = filter._lower_level;
	this->_lower_interval = filter._lower_interval;
	this->_upper_level = filter._upper_level;
	this->_ratio = filter._ratio;

	this->_mer_map = filter._mer_map;
	this->_debug = filter._debug;
}

Filter::Filter() {
	this->_mer_length = 0;
	this->_lower_level = 0;
	this->_lower_interval = 0;
	this->_upper_level = 0;
	this->_ratio = 0.;
	this->_debug = false;
}

bool Filter::insertMer(const Read& read, int score) throw(MerLengthError) {
	if (this->_mer_length == 0) {
		this->_mer_length = read.size();
	} else if (read.size() != this->_mer_length) {
		throw MerLengthError();
	} else if (score <= this->_lower_level) {
		// score <= low level
		return false;
	}

	this->_mer_map.insert(mer_map::value_type(read, score));
	return true;
}

bool Filter::insertMers(Fasta& fasta) {
	bool retval(false);
	while (!fasta.eof()) {
		const Fasta::Item item(fasta.getItem());
		int score(0);
		std::string str;
		try {
			str = item.getInfo();
			score = boost::lexical_cast<int>(str);
			const Read read(item.getRead());
			const bool flg(insertMer(read, score));
			retval = retval || flg;
		} catch(const boost::bad_lexical_cast& e) {
			if (_debug) {
				std::cerr << e.what() << std::endl;
			}
			continue;
		} catch(const MerLengthError& e) {
			continue;
		}
	}
	return retval;
}

bool Filter::join(const Filter& filter) throw(MerLengthError, LowLevelError) {
	if (filter._lower_level != this->_lower_level) {
		throw LowLevelError();
	}
	if (this->_mer_length == 0) {
		this->_mer_length = filter._mer_length;
	} else if (filter._mer_length != this->_mer_length) {
		throw MerLengthError();
	}

	this->_mer_map.insert(filter._mer_map.begin(), filter._mer_map.end());
	return true;
}

std::vector<unsigned int> Filter::scores(const Read& read) const {
	std::vector<unsigned int> retval;
	const int length(read.size() - _mer_length + 1);
	if (_mer_length == 0 || length <= 0) {
		return retval;
	}
	for (int i(0); i < length; i++) {
		const Read sub(read.sub(i, _mer_length));
		if (!sub.isDefinite()) {
			continue;
		}
		retval.push_back(_getScore(sub,0));
	}
	return retval;
}

bool Filter::check(std::vector<unsigned int> scores) const {
	if (scores.size() == 0) {
		return false;
	}

	int lower_total(0), lower_count(0);
	int upper_total(0), upper_count(0);

	for (std::vector<unsigned int>::const_iterator itr(scores.begin());
			itr != scores.end(); itr++) {
		if (*itr <= _lower_level) {
			lower_count++;
			lower_total += *itr;
		} else {
			upper_count++;
			upper_total += *itr;
			if (lower_count < _lower_interval) {
				lower_count = 0;
			}
		}
	}

	if (upper_count == 0) {
		return true;
	}
	if (lower_count < _lower_interval) {
		return true;
	}

	double upper_average(double(upper_total)/upper_count),
		   lower_average(double(lower_total)/lower_count);
	return upper_average < _upper_level || upper_average < lower_average * _ratio;
}

bool Filter::check(const Read& read) const {
	return check(scores(read));
}

int Filter::average(const Read& read) const {
	std::vector<unsigned int> scores(this->scores(read));
	if (scores.size() == 0) {
		return 0;
	}

	unsigned int total(0);
	for (std::vector<unsigned int>::const_iterator itr(scores.begin());
			itr != scores.end(); itr++) {
		total += *itr;
	}

	return total/scores.size();
}

int Filter::_getScore(const Read& read, const int default_value) const
		throw(MerLengthError){
	if (read.size() != _mer_length) {
		throw MerLengthError();
	}
	int score(0);
	mer_map::const_iterator itr(_mer_map.find(read));
	if (itr != _mer_map.end()) {
		score = (*itr).second;
	} else {
		mer_map::const_iterator comp(_mer_map.find(read.complement()));
		if (comp != _mer_map.end()) {
			score = (*comp).second;
		} else {
			return default_value;
		}
	}
	return score;
}

void Filter::setDebug(bool debug) {
	_debug = debug;
}

