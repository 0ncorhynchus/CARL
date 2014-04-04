// filter.hpp
// written by S.Kato

#ifndef __FILTER_hpp
#define __FILTER_hpp

#include <string>
#include <stdexcept>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <unordered_map>
#include "read.hpp"
#include "fasta.hpp"

namespace carl {

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

    typedef unsigned int score_type;
    typedef std::unordered_map<Read, score_type> map_type;

private:
    map_type _mer_map;
    Read::size_type  _mer_length;
    score_type _lower_level;
    unsigned int _lower_interval;
    double _ratio;
    int _getScore(const Read& read, const score_type default_value) const
        throw(MerLengthError);

public:
    Filter(score_type lower_level, unsigned int lower_interval, double ratio);
    Filter(const Filter& filter);
    Filter();
    bool insertMer(const Read& read, score_type score) throw(MerLengthError);
    bool insertMers(Fasta& fasta);
    bool join(const Filter& filter) throw(MerLengthError, LowerLevelError);
    std::vector<score_type> scores(const Read& read) const;
    bool check(std::vector<score_type> scores) const;
    bool check(const Read& read) const;
    double average(std::vector<score_type> scores) const;
    double average(const Read& read) const;
    int size() const {
        return this->_mer_map.size();
    }
};

} // carl

#endif
