// filter.hpp
// written by S.Kato

#ifndef __FILTER_hpp
#define __FILTER_hpp

#include <mutex>
#include <string>
#include <map>
#include "fasta.hpp"

typedef std::map<Read, int> mer_map;
typedef mer_map::iterator mer_iterator;

class Filter {
private:
	mutable std::mutex m_mutex;
	mer_map _mer_map;
	int _mer_length;
	int _low_level, _top_level;
	int _low_interval;
	bool _debug;
	int _getScore(Read read);
public:
	Filter(int mer_length);
	Filter(int mer_length, int low_level, int low_interval, int top_level);
	void setMerLength(int length);
	void addMer(Read read, int score);
	bool check(int scores[]);
	mer_iterator begin();
	mer_iterator end();
	void getScoreList(Read read, int scores[]);
	void setDebug(bool);
};

#endif
