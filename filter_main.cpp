#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "filter.hpp"

void insertMer(const std::string mer_file, Filter& filter) {
    Fasta mers(mer_file);
    filter.insertMers(mers);
}

int average(const std::string infile, std::ostream& str, const Filter& filter) {
    Fasta fasta(infile);
    while (!fasta.eof()) {
        const std::pair<std::string, std::string> item(fasta.getItemStrings());
        const Read read(item.second);

        if (read.size() == 0)
            continue;

        const std::string info(item.first);
        const double average(filter.average(read));
        str << ">" << info << std::endl;
        str << average << std::endl;
    }
    int score(0);
    return score;
}

void check(const std::string infile, std::ostream& str, const Filter& filter) {
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

void filter(const std::string& read_file, const std::string& mers_file,
        const unsigned int& lower_level, const unsigned int& low_interval,
        const double& ratio, const unsigned int& cpua, const unsigned int& cpub) {

    Filter *filter = new Filter(lower_level, low_interval, ratio);

    /*
     * Importing mer from a file
     */
    boost::uuids::random_generator rng;
    const boost::uuids::uuid id = rng();
    std::string identifier(boost::lexical_cast<std::string>(id));

    if (cpua == 1) {
        insertMer(mers_file, *filter);
    } else {
        std::string *filenames = new std::string[cpua];
        Filter *filter_set = new Filter[cpua];
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
        for (int i(0); i < cpua; i++)
        {
            std::ostringstream oss;
            oss << "/tmp/filter_read_splita_" << identifier << "_" << i;
            filenames[i] = oss.str();

            std::ofstream ofs(filenames[i]);
            for (int j(lines*i/cpua); j < lines*(i+1)/cpua; j++) {
                getline(ifs, buf);
                ofs << buf << std::endl;
                getline(ifs, buf);
                ofs << buf << std::endl;
            }
            ofs.close();

            filter_set[i] = Filter(lower_level, low_interval, ratio);
            threads.create_thread(boost::bind(&insertMer, filenames[i] , std::ref(filter_set[i])));
        }
        ifs.close();
        threads.join_all();

        for (int i(0); i < cpua; i++) {
            try {
                filter->join(filter_set[i]);
            } catch(const Filter::LowerLevelError& e) {
                std::cerr << e.what() << std::endl;
            } catch(const Filter::MerLengthError& e) {
                std::cerr << e.what() << std::endl;
            }
            remove(filenames[i].c_str());
        }
        delete[] filenames;
        delete[] filter_set;
    }

    if (cpub == 1) {
        check(read_file, std::cout, *filter);
    } else {
        std::string *infiles = new std::string[cpub];
        std::string *outfiles = new std::string[cpub];
        std::ofstream *outstreams = new std::ofstream[cpub];
        Filter *filter_set = new Filter[cpub];

        boost::thread_group threads;
        std::ifstream ifs(read_file);
        std::string buf;
        int lines(0);
        while (ifs && getline(ifs, buf)) {
            lines++;
        }
        lines /= 2;
        ifs.close();

        ifs.open(read_file);
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
                        std::ref(outstreams[i]), std::ref(filter_set[i])));
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

void calculate_average(const std::string& read_file, const std::string& mers_file,
        const unsigned int& cpua, const unsigned int& cpub) {

    Filter *filter = new Filter(1,0,0);

    /*
     * Importing mer from a file
     */
    boost::uuids::random_generator rng;
    const boost::uuids::uuid id = rng();
    std::string identifier(boost::lexical_cast<std::string>(id));

    if (cpua == 1) {
        insertMer(mers_file, *filter);
    } else {
        std::string *filenames = new std::string[cpua];
        Filter *filter_set = new Filter[cpua];
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
        for (int i(0); i < cpua; i++)
        {
            std::ostringstream oss;
            oss << "/tmp/filter_read_splita_" << identifier << "_" << i;
            filenames[i] = oss.str();

            std::ofstream ofs(filenames[i]);
            for (int j(lines*i/cpua); j < lines*(i+1)/cpua; j++) {
                getline(ifs, buf);
                ofs << buf << std::endl;
                getline(ifs, buf);
                ofs << buf << std::endl;
            }
            ofs.close();

            filter_set[i] = Filter(1,0,0);
            threads.create_thread(boost::bind(&insertMer, filenames[i] , std::ref(filter_set[i])));
        }
        ifs.close();
        threads.join_all();

        for (int i(0); i < cpua; i++) {
            try {
                filter->join(filter_set[i]);
            } catch(const Filter::LowerLevelError& e) {
                std::cerr << e.what() << std::endl;
            } catch(const Filter::MerLengthError& e) {
                std::cerr << e.what() << std::endl;
            }
            remove(filenames[i].c_str());
        }
        delete[] filenames;
        delete[] filter_set;
    }

    if (cpub == 1) {
        average(read_file, std::cout, *filter);
    } else {
        std::string *infiles = new std::string[cpub];
        std::string *outfiles = new std::string[cpub];
        std::ofstream *outstreams = new std::ofstream[cpub];
        Filter *filter_set = new Filter[cpub];

        boost::thread_group threads;
        std::ifstream ifs(read_file);
        std::string buf;
        int lines(0);
        while (ifs && getline(ifs, buf)) {
            lines++;
        }
        lines /= 2;
        ifs.close();

        ifs.open(read_file);
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
            threads.create_thread(boost::bind(&average, infiles[i],
                        std::ref(outstreams[i]), std::ref(filter_set[i])));
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

int main(int argc, char** argv) {
    std::string command(argv[0]);
    std::string usage("usage: " + command + " read_file mer_file [options]");

    unsigned int lower_level(0), low_interval(0);
    double ratio;
    int result;
    unsigned int cpua(1), cpub(1);
    using namespace boost::program_options;
    options_description options0(""), options1("");
    options0.add_options()
        (",f", value<unsigned int>(&lower_level), "lower_level")
        (",m", value<unsigned int>(&low_interval), "low_frequence")
        (",r", value<double>(&ratio), "ratio")
        (",a", value<unsigned int>(&cpua)->default_value(1), "threads for creating maps")
        (",b", value<unsigned int>(&cpub)->default_value(1), "threads for calculating");
    options1.add_options()
        ("average", "calculating average scores");
    options0.add(options1);

    if (argc < 3) {
        std::cerr << usage << std::endl;
        std::cerr << "OPTIONS" << std::endl;
        std::cerr << options0 << std::endl;
        return 1;
    }
    std::string read_file(argv[1]);
    std::string mers_file(argv[2]);

    variables_map values;
    try {
        store(parse_command_line(argc, argv, options0), values);
        notify(values);
        if (values.count("average")) {
            calculate_average(read_file, mers_file, cpua, cpub);
        } else {
            filter(read_file, mers_file, lower_level, low_interval, ratio, cpua, cpub);
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
    }
}
