#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include "filter.hpp"

void insertMer(const std::string filename, Filter& filter, bool debug) {
	Fasta mers(filename);
	while (!mers.eof()) {
		FastaItem item(mers.getItem());
		int score = 0;
		std::string str;
		try {
			str = item.getInfo();
			score = boost::lexical_cast<int>(str);
		} catch(boost::bad_lexical_cast) {
			if (debug) {
				std::cerr << "boost::bad_lexical_cast" << str << std::endl;
			}
			continue;
		}
		Read read = item.getRead();
		bool flg = filter.insertMer(read, score);
		if (!flg && debug) {
			std::cerr << "Not Import mer : " << read.tostring() << " @" << score << std::endl;
		}
	}
}

void check(const std::string infile, std::ostream& str, const Filter& filter, bool debug) {
	std::cerr << &filter << std::endl;
	Fasta fasta(infile);
	while (!fasta.eof()) {
		std::pair<std::string, std::string> item(fasta.getItemStrings());
		Read read(item.second);

		if (read.size() == 0)
			continue;

		if (filter.check(read)) {
			std::string info = item.first;
			std::string seq = item.second;
			str << ">" << info << std::endl;
			str << seq << std::endl;
		}
	}
}

int main(int argc, char** argv) {
	std::string command(argv[0]);
	std::string usage(
			"usage: " + command + " filename mer_file" +
			" [-d] [-f low_level] [-m low_frequence] [-t top_level] [-a threads] [-b threads]"
			);
	if (argc < 3) {
		std::cerr << usage << std::endl;
		return 1;
	}
	std::string filename(argv[1]);
	std::string mers_file(argv[2]);

	unsigned int top_level, low_level(0), low_interval;
	int result;
	unsigned int cpus(1), cpub(1);
	bool debug(false);
	while ((result = getopt(argc, argv, "df:m:t:a:b:")) != -1) {
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
				case 'b':
					cpub = boost::lexical_cast<unsigned int>(optarg);
					if (cpub == 0) {
						cpub = 1;
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

	std::ostringstream oss;
	int mers_sindex(mers_file.find_last_of("/") + 1);
	int mers_eindex(mers_file.find_last_of("."));
	int mers_sublen(0);
	if (mers_eindex == std::string::npos || mers_eindex < mers_sindex) {
		mers_sublen = mers_file.length() - mers_sindex;
	} else {
		mers_sublen = mers_eindex - mers_sindex;
	}
	int file_sindex(filename.find_last_of("/") + 1);
	int file_eindex(filename.find_last_of("."));
	int file_sublen(0);
	if (file_eindex == std::string::npos || file_eindex < file_sindex) {
		file_sublen = filename.length() - file_sindex;
	} else {
		file_sublen = file_eindex - file_sindex;
	}
	oss << mers_file.substr(mers_sindex, mers_sublen) << "_" << filename.substr(file_sindex, file_sublen);
	oss <<  "_" << low_level << "_" << low_interval << "_" << top_level << std::flush;
	std::string identifier(oss.str());


	if (cpus == 1) {
		insertMer(mers_file, *filter, debug);
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
		lines /= 2;
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
			oss << "/tmp/filter_read_splita_" << identifier << "_" << i;
			filenames[i] = oss.str();

			std::ofstream ofs(filenames[i]);
			for (int j(lines*i/cpus); j < lines*(i+1)/cpus; j++) {
				getline(ifs, buf);
				ofs << buf << std::endl;
				getline(ifs, buf);
				ofs << buf << std::endl;
			}
			ofs.close();

			filter_set[i] = Filter(low_level, low_interval, top_level);
			filter_set[i].setDebug(debug);
			threads.create_thread(boost::bind(&insertMer, filenames[i] , std::ref(filter_set[i]), debug));
		}
		ifs.close();
		threads.join_all();

		bool isJoined(true);
		for (int i(0); i < cpus; i++) {
			isJoined = isJoined && filter->insertMers(filter_set[i]);
			remove(filenames[i].c_str());
		}
		if (!isJoined && debug) {
			std::cerr << "Failed join maps" << std::endl;
		}
		delete[] filenames;
		delete[] filter_set;
	}

	if (debug) {
		std::cerr << "finish map" << std::endl;
		std::cerr << "map size: " <<  filter->size() << std::endl;
	}

	if (cpub == 1) {
		check(filename, std::cout, *filter, debug);
	} else {
		std::string *infiles = new std::string[cpub];
		std::string *outfiles = new std::string[cpub];
		std::ofstream *outstreams = new std::ofstream[cpub];
		Filter *filter_set = new Filter[cpub];

		boost::thread_group threads;
		std::ifstream ifs(filename);
		std::string buf;
		int lines(0);
		while (ifs && getline(ifs, buf)) {
			lines++;
		}
		lines /= 2;
		ifs.close();

		ifs.open(filename);
		for (int i(0); i < cpub; i++)
		{
			std::ostringstream ioss, ooss;
			ioss << "/tmp/filter_read_splitb_" << identifier <<  "_" << i;
			ooss << "/tmp/filter_read_split_out_" << identifier << "_" << i;
			infiles[i] = ioss.str();
			outfiles[i] = ooss.str();

			std::ofstream ofs(infiles[i]);
			for (int j(lines*i/cpub); j < lines*(i+1)/cpub; j++) {
				getline(ifs, buf);
				ofs << buf << std::endl;
				getline(ifs, buf);
				ofs << buf << std::endl;
			}
			ofs.close();

			outstreams[i].open(outfiles[i].c_str());
			filter_set[i] = Filter(*filter);
			threads.create_thread(boost::bind(&check, infiles[i],
						std::ref(outstreams[i]), std::ref(filter_set[i]), debug));
		}
		ifs.close();
		threads.join_all();

		for (int i(0); i < cpub; i++) {
			outstreams[i].close();
			std::ifstream tmpifs(outfiles[i].c_str());
			std::string buffer;
			while (tmpifs && !tmpifs.eof()) {
				if (getline(tmpifs, buffer)) {
					std::cout << buffer << std::endl;
				}
			}
			tmpifs.close();
			remove(infiles[i].c_str());
			remove(outfiles[i].c_str());
		}
		delete[] filter_set;
		delete[] outstreams;
		delete[] infiles;
		delete[] outfiles;
	}

	delete filter;
}
