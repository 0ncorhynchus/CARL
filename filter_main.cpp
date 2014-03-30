#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "filter.hpp"

void import_mer(const std::string mer_file, Filter& filter) {
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

Filter import_mer_with_multi_thread(const std::string mers_file,
        const Filter& parent, const int num_thread, const std::string identifier) {
    Filter retval(parent);
    if (num_thread <= 1) {
        import_mer(mers_file, retval);
    } else {
        std::vector<std::string> filenames;
        std::vector<Filter> filters;
        filenames.reserve(num_thread);
        filters.reserve(num_thread);

        boost::thread_group threads;
        std::ifstream ifs(mers_file);
        std::string buff;
        int lines(0);
        while (ifs && getline(ifs, buff)) {
            lines++;
        }
        lines /= 2;
        ifs.close();

        ifs.open(mers_file);
        for (int i(0); i < num_thread; i++) {
            std::ostringstream oss;
            oss << "/tmp/filter_mer_split_" << identifier << "_" << i;
            filenames.push_back(oss.str());

            std::ofstream ofs(filenames.back());
            for (int j(lines*i/num_thread); j < lines*(i+1)/num_thread; j++) {
                getline(ifs, buff);
                ofs << buff << std::endl;
                getline(ifs, buff);
                ofs << buff << std::endl;
            }
            ofs.close();

            filters.push_back(parent);
            threads.create_thread(boost::bind(&import_mer, filenames.back(),
                        std::ref(filters.back())));
        }
        ifs.close();
        threads.join_all();

        for (int i(0); i < num_thread; i++) {
            try {
                retval.join(filters.at(i));
                remove(filenames.at(i).c_str());
            } catch(const Filter::LowerLevelError& e) {
                std::cerr << e.what() << std::endl;
                std::ostringstream oss;
                oss << "./not_imported_" << identifier << "_" << i << ".fasta";
                rename(filenames.at(i).c_str(), oss.str().c_str());
            } catch(const Filter::MerLengthError& e) {
                std::cerr << e.what() << std::endl;
                std::ostringstream oss;
                oss << "./not_imported_" << identifier << "_" << i << ".fasta";
                rename(filenames.at(i).c_str(), oss.str().c_str());
            }
        }
    }
    return retval;
}

void filter(const std::string& read_file, const std::string& mers_file,
        const unsigned int& lower_level, const unsigned int& low_interval,
        const double& ratio, const unsigned int& cpua, const unsigned int& cpub) {

    Filter filter(lower_level, low_interval, ratio);

    /*
     * Importing mer from a file
     */
    boost::uuids::random_generator rng;
    const boost::uuids::uuid id = rng();
    const std::string identifier(boost::lexical_cast<std::string>(id));

    filter = import_mer_with_multi_thread(mers_file, filter, cpua, identifier);

    if (cpub == 1) {
        check(read_file, std::cout, filter);
    } else {
        std::vector<std::string> infiles;
        std::vector<std::string> outfiles;
        std::ofstream* ostreams = new std::ofstream[cpub];
        std::vector<Filter> filters;
        infiles.reserve(cpub);
        outfiles.reserve(cpub);
        filters.reserve(cpub);

        boost::thread_group threads;
        std::ifstream ifs(read_file);
        std::string buff;
        int lines(0);
        while (ifs && getline(ifs, buff)) {
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
            infiles.push_back(ioss.str());
            outfiles.push_back(ooss.str());

            std::ofstream ofs(infiles.back());
            for (int j(lines*i/cpub); j < lines*(i+1)/cpub; j++) {
                getline(ifs, buff);
                ofs << buff << std::endl;
                getline(ifs, buff);
                ofs << buff << std::endl;
            }
            ofs.close();

            ostreams[i].open(outfiles.back().c_str());
            filters.push_back(filter);
            threads.create_thread(boost::bind(&check, infiles.back(),
                        std::ref(ostreams[i]), std::ref(filters.back())));
        }
        ifs.close();
        threads.join_all();

        for (int i(0); i < cpub; i++) {
            ostreams[i].close();
            std::ifstream tmpifs(outfiles.at(i).c_str());
            std::string buffer;
            while (tmpifs && !tmpifs.eof()) {
                if (getline(tmpifs, buffer)) {
                    std::cout << buffer << std::endl;
                }
            }
            tmpifs.close();
            remove(infiles.at(i).c_str());
            remove(outfiles.at(i).c_str());
        }
        delete[] ostreams;
    }

}

void calculate_average(const std::string& read_file, const std::string& mers_file,
        const unsigned int& cpua, const unsigned int& cpub) {

    Filter filter(1,0,0);

    /*
     * Importing mer from a file
     */
    boost::uuids::random_generator rng;
    const boost::uuids::uuid id = rng();
    const std::string identifier(boost::lexical_cast<std::string>(id));

    filter = import_mer_with_multi_thread(mers_file, filter, cpua, identifier);

    if (cpub == 1) {
        average(read_file, std::cout, filter);
    } else {
        std::vector<std::string> infiles, outfiles;
        std::ofstream* ostreams = new std::ofstream[cpub];
        std::vector<Filter> filters;
        infiles.reserve(cpub);
        outfiles.reserve(cpub);
        filters.reserve(cpub);

        boost::thread_group threads;
        std::ifstream ifs(read_file);
        std::string buff;
        int lines(0);
        while (ifs && getline(ifs, buff)) {
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
            infiles.push_back(ioss.str());
            outfiles.push_back(ooss.str());

            std::ofstream ofs(infiles.back());
            for (int j(lines*i/cpub); j < lines*(i+1)/cpub; j++) {
                getline(ifs, buff);
                ofs << buff << std::endl;
                getline(ifs, buff);
                ofs << buff << std::endl;
            }
            ofs.close();

            ostreams[i].open(outfiles.back().c_str());
            filters.push_back(filter);
            threads.create_thread(boost::bind(&average, infiles.back(),
                        std::ref(ostreams[i]), std::ref(filters.back())));
        }
        ifs.close();
        threads.join_all();

        for (int i(0); i < cpub; i++) {
            ostreams[i].close();
            std::ifstream tmpifs(outfiles.at(i).c_str());
            std::string buffer;
            while (tmpifs && !tmpifs.eof()) {
                if (getline(tmpifs, buffer)) {
                    std::cout << buffer << std::endl;
                }
            }
            tmpifs.close();
            remove(infiles.at(i).c_str());
            remove(outfiles.at(i).c_str());
        }
        delete[] ostreams;
    }
}

int main(int argc, char** argv) {
    std::string command(argv[0]);
    std::string usage("usage: " + command + " read_file mer_file [options]");

    unsigned int lower_level(0), low_interval(0);
    double ratio;
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
