// filter.hpp
// written by S.Kato

#ifndef __FILTER_hpp
#define __FILTER_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>
#include "fasta.hpp"

typedef boost::unordered_map<Read, int> mer_map;

class Filter {
public:
	class MerLengthError : public std::domain_error {
	public:
		MerLengthError() :
			std::domain_error::domain_error("MerLengthError")
		{
		}
	};

	class LowLevelError : public std::domain_error {
	public:
		LowLevelError() :
			std::domain_error::domain_error("LowLevelError")
		{
		}
	};

private:
	mer_map _mer_map;
	int _mer_length;
	int _lower_level, _upper_level;
	int _lower_interval;
	double _ratio;
	bool _debug;
	int _getScore(const Read& read, const int default_value) const
		throw(MerLengthError);

public:
	Filter(int lower_level, int lower_interval, int upper_level, double ratio);
	Filter(const Filter& filter);
	Filter();
	bool insertMer(const Read& read, int score) throw(MerLengthError);
	bool insertMers(Fasta& fasta);
	bool join(const Filter& filter) throw(MerLengthError, LowLevelError);
	std::vector<unsigned int> scores(const Read& read) const;
	bool check(std::vector<unsigned int> scores) const;
	bool check(const Read& read) const;
	int average(const Read& read) const;
	void setDebug(bool);
	int size() {
		return this->_mer_map.size();
	}
};

#endif
