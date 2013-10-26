// filter.cpp
// written by S.Kato

#include <iostream>
#include <string>
#include "filter.hpp"

Filter::Filter(int low_level, int low_interval, int top_level) {
	this->_mer_length = 0;
	this->_low_level = low_level;
	this->_low_interval = low_interval;
	this->_top_level = top_level;
}

Filter::Filter(const Filter& filter) {
	this->_mer_length = filter._mer_length;
	this->_low_level = filter._low_level;
	this->_low_interval = filter._low_interval;
	this->_top_level = filter._top_level;
	this->_mer_map = filter._mer_map;
	this->_debug = filter._debug;
}

Filter::Filter() {
	this->_mer_length = 0;
	this->_low_level = 0;
	this->_low_interval = 0;
	this->_top_level = 0;
}

bool Filter::insertMer(const Read read, int score) {
	if (this->_mer_length == 0) {
		this->_mer_length = read.size();
	} else if (read.size() != this->_mer_length || score <= this->_low_level) {
		return false;
	}

	this->_mer_map.insert(mer_map::value_type(read, score));
	return true;
}

bool Filter::insertMers(const Filter& filter) {
	if (filter._low_level != this->_low_level)
		return false;
	if (this->_mer_length == 0) {
		this->_mer_length = filter._mer_length;
	} else if (filter._mer_length != this->_mer_length) {
		return false;
	}

	this->_mer_map.insert(filter._mer_map.begin(), filter._mer_map.end());
	return true;
}

bool Filter::check(Read read) const {
	if (this->_mer_length == 0 || read.size() < this->_mer_length) {
		return false;
	}
	int low_count(0), total(0), count(0);
	for (int i(0); i < read.size() - this->_mer_length; i++) {
		Read sub(read.sub(i, this->_mer_length));
		if (!read.isDefinite()) {
			continue;
		}
		int score(this->_getScore(sub));
		if (score <= this->_low_level) {
			low_count++;
		} else {
			count++;
			total += score;
			if (low_count < this->_low_interval) {
				low_count = 0;
			}
		}
	}

	if (count == 0) {
		if (_debug) {
			std::cerr << "# low_count=" << low_count << ", count = 0" << std::endl;
		}
		return false;
	}

	if (_debug) {
		std::cerr << "# low_count=" << low_count << ", total=" << total << std::endl;
	}
	return low_count < this->_low_interval || total < _top_level * count;
}

int Filter::_getScore(Read read) const {
	int score(0);
	mer_map::const_iterator itr(_mer_map.find(read));
	if (itr != _mer_map.end()) {
		score = (*itr).second;
	} else {
		mer_map::const_iterator comp(_mer_map.find(read.complement()));
		if (comp != _mer_map.end()) {
			score = (*comp).second;
		}
	}
	return score;
}

void Filter::setDebug(bool debug) {
	_debug = debug;
}

