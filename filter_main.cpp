#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include "filter.hpp"

void run(const std::string filename, Filter& filter) {
	Fasta mers(filename);
	while (!mers.eof()) {
		FastaItem item(mers.getItem());
		int score = 0;
		try {
			std::string str(item.getInfo());
			score = boost::lexical_cast<int>(str);
		} catch(boost::bad_lexical_cast) {
			continue;
		}
		filter.insertMer(item.getRead(), score);
	}
}

int main(int argc, char** argv) {
	std::string command(argv[0]);
	std::string usage(
			"usage: " + command + " filename mer_file" +
			" [-d] [-f low_level] [-m low_frequence] [-t top_level] [-a threads]"
			);
	if (argc < 3) {
		std::cerr << usage << std::endl;
		return 1;
	}
	std::string filename(argv[1]);
	std::string mers_file(argv[2]);

	unsigned int top_level, low_level, low_interval;
	int result;
	unsigned int cpus(1);
	bool debug(false);
	while ((result = getopt(argc, argv, "df:m:t:a:")) != -1) {
		try {
			switch(result) {
				case 'd':
					debug = true;
					break;
				case 'f':
					low_level = boost::lexical_cast<unsigned int>(optarg);
					break;
				case 'm':
					low_interval = boost::lexical_cast<unsigned int>(optarg);
					break;
				case 't':
					top_level = boost::lexical_cast<unsigned int>(optarg);
					break;
				case 'a':
					cpus = boost::lexical_cast<unsigned int>(optarg);
					if (cpus == 0) {
						cpus = 1;
					}
					break;
				case ':':
				case '?':
					std::cerr << usage << std::endl;
					return 1;
			}
		} catch(boost::bad_lexical_cast) {
			std::cerr << usage << std::endl;
			return 1;
		}
	}

	Filter *filter = new Filter(low_level, low_interval, top_level);
	filter->setDebug(debug);

	/*
	 * Importing mer from a file
	 */
	if (cpus == 1) {
		run(mers_file, *filter);
	} else {
		std::string *filenames = new std::string[cpus];
		Filter *filter_set = new Filter[cpus];
		boost::thread_group threads;
		std::ifstream ifs(mers_file);
		std::string buf;
		int lines(0);
		while (ifs && getline(ifs, buf)) {
			lines++;
		}
		ifs.close();

		ifs.open(mers_file);
		for (int i(0); i < cpus; i++)
		{
			std::ostringstream oss;
			int sindex(mers_file.find_last_of("/") + 1);
			int eindex(mers_file.find_last_of("."));
			int sublen(0);
			if (eindex == std::string::npos || eindex < sindex) {
				sublen = mers_file.length() - sindex;
			} else {
				sublen = eindex - sindex;
			}
			oss << "/tmp/filter_read_split_" << mers_file.substr(sindex, sublen) << "_" << low_level << "_" << low_interval << "_" << top_level << "_" << i;
			filenames[i] = oss.str();

			std::ofstream ofs(filenames[i]);
			for (int j(lines*i/cpus); j < lines*(i+1)/cpus; j++) {
				getline(ifs, buf);
				ofs << buf << std::endl;
			}
			ofs.close();

			filter_set[i] = Filter(low_level, low_interval, top_level);
			threads.create_thread(boost::bind(&run, filenames[i] , std::ref(filter_set[i])));
		}
		ifs.close();
		threads.join_all();

		for (int i(0); i < cpus; i++) {
			filter->insertMers(filter_set[i]);
			remove(filenames[i].c_str());
		}
		delete[] filenames;
		delete[] filter_set;
	}

	if (debug)
		std::cerr << "finish map" << std::endl;

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
