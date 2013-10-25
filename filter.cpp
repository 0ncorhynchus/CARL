// filter.cpp
// written by S.Kato

#include <iostream>
#include <string>

#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include "filter.hpp"

Filter::Filter(int low_level, int low_interval, int top_level) {
	this->_mer_length = 0;
	this->_low_level = low_level;
	this->_low_interval = low_interval;
	this->_top_level = top_level;
}

bool Filter::addMer(const Read read, int score) {
	if (this->_mer_length == 0) {
		this->_mer_length = read.size();
	} else if (read.size() != this->_mer_length || score <= this->_low_level) {
		return false;
	}

	this->_mer_map[read] = score;
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
			std::cout << "# low_count=" << low_count << ", count = 0" << std::endl;
		}
		return false;
	}

	if (_debug) {
		std::cout << "# low_count=" << low_count << ", total=" << total << std::endl;
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

int main(int argc, char** argv) {
	std::string command(argv[0]);
	std::string usage(
			"usage: " + command + " filename mer_file" +
			" [-d] [-f low_level] [-m low_frequence] [-t top_level] [-a threads]"
			);
	if (argc < 3) {
		std::cout << usage << std::endl;
		return 1;
	}
	std::string filename(argv[1]);
	std::string mers_file(argv[2]);

	int top_level, low_level, low_interval;
	int result;
	int cpus(1);
	bool debug(false);
	while ((result = getopt(argc, argv, "df:m:t:a:")) != -1) {
		try {
			switch(result) {
				case 'd':
					debug = true;
					break;
				case 'f':
					low_level = boost::lexical_cast<int>(optarg);
					break;
				case 'm':
					low_interval = boost::lexical_cast<int>(optarg);
					break;
				case 't':
					top_level = boost::lexical_cast<int>(optarg);
					break;
				case 'a':
					cpus = boost::lexical_cast<int>(optarg);
					break;
				case ':':
				case '?':
					std::cout << usage << std::endl;
					return 1;
			}
		} catch(boost::bad_lexical_cast) {
			std::cout << usage << std::endl;
			return 1;
		}
	}

	Fasta *mers = new Fasta(mers_file);
	Filter *filter = new Filter(low_level, low_interval, top_level);
	filter->setDebug(debug);

	/*
	 * Importing mer from a file
	 */
	/*
	boost::thread_group threads;
	for (int i(0); i < cpus; i++)
	{
		threads.create_thread([&] {
					while (!mers->eof()) {
						FastaItem item = mers->getItem();
						int score = 0;
						try {
							std::string str = item.getInfo();
							score = boost::lexical_cast<int>(str);
						} catch(boost::bad_lexical_cast) {
							continue;
						}

						filter->addMer(item.getRead(), score);
					}
				});
	}
	threads.join_all();
	*/

	while (!mers->eof()) {
		FastaItem item(mers->getItem());
		int score = 0;
		try {
			std::string str(item.getInfo());
			score = boost::lexical_cast<int>(str);
		} catch(boost::bad_lexical_cast) {
			continue;
		}
		filter->addMer(item.getRead(), score);
	}
	delete mers;

	if (debug)
		std::cout << "finish map" << std::endl;

	Fasta *fasta = new Fasta(filename);
	while (! fasta->eof()) {
		FastaItem item(fasta->getItem());
		Read read(item.getRead());

		if (read.size() == 0)
			continue;

		if (filter->check(read)) {
			std::string info = item.getInfo();
			std::cout << ">" << info << std::endl;
			std::cout << read.tostring() << std::endl;
		}
	}

	delete fasta;
	delete filter;
}
