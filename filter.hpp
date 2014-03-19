// filter.hpp
// written by S.Kato

#ifndef __FILTER_hpp
#define __FILTER_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <unordered_map>
#include "fasta.hpp"

typedef std::unordered_map<Read, int> mer_map;

class Filter {
public:
    class MerLengthError : public std::domain_error {
    public:
        MerLengthError(const std::string& what_arg) :
            std::domain_error::domain_error("MerLengthError: " + what_arg)
        {
        }
    };

    class LowerLevelError : public std::domain_error {
    public:
        LowerLevelError(const std::string& what_arg) :
            std::domain_error::domain_error("LowerLevelError: " + what_arg)
        {
        }
    };

private:
    mer_map _mer_map;
    int _mer_length;
    int _lower_level;
    int _lower_interval;
    double _ratio;
    int _getScore(const Read& read, const int default_value) const
        throw(MerLengthError);

public:
    Filter(int lower_level, int lower_interval, double ratio);
    Filter(const Filter& filter);
    Filter();
    bool insertMer(const Read& read, int score) throw(MerLengthError);
    bool insertMers(Fasta& fasta);
    bool join(const Filter& filter) throw(MerLengthError, LowerLevelError);
    std::vector<unsigned int> scores(const Read& read) const;
    bool check(std::vector<unsigned int> scores) const;
    bool check(const Read& read) const;
    double average(std::vector<unsigned int> scores) const;
    double average(const Read& read) const;
    int size() const {
        return this->_mer_map.size();
    }
};

#endif
