// filter.cpp
// written by S.Kato

#include <iostream>
#include <string>
#include "filter.hpp"

namespace carl {

Filter::Filter(score_type lower_level, unsigned int lower_interval, double ratio) {
    this->_mer_length = 0;
    this->_lower_level = lower_level;
    this->_lower_interval = lower_interval;
    this->_ratio = ratio;
}

Filter::Filter(const Filter& filter) {
    this->_mer_length = filter._mer_length;
    this->_lower_level = filter._lower_level;
    this->_lower_interval = filter._lower_interval;
    this->_ratio = filter._ratio;

    this->_mer_map = filter._mer_map;
}

Filter::Filter() {
    this->_mer_length = 0;
    this->_lower_level = 0;
    this->_lower_interval = 0;
    this->_ratio = 0.;
}

bool Filter::insertMer(const Read& read, score_type score) throw(MerLengthError) {
    if (!read.isDefinite())
        return false;
    if (this->_mer_length == 0) {
        this->_mer_length = read.size();
    } else if (read.size() != this->_mer_length) {
        std::ostringstream oss;
        oss << _mer_length << " is not " << read.size();
        oss << ", Failed inserting " << read.tostring();
        throw MerLengthError(oss.str());
    } else if (score <= this->_lower_level) {
        return false;
    }

    this->_mer_map.insert(map_type::value_type(read, score));
    return true;
}

bool Filter::insertMers(Fasta& fasta) {
    bool retval(false);
    while (!fasta.eof()) {
        const Fasta::Item item(fasta.getItem());
        int score(0);
        std::string str;
        try {
            str = item.getInfo();
            score = boost::lexical_cast<int>(str);
            const Read read(item.getRead());
            const bool flg(insertMer(read, score));
            retval = retval || flg;
        } catch(const boost::bad_lexical_cast& e) {
            std::cerr << e.what() << ", from \"" << str << "\" to <int>" << std::endl;
            continue;
        } catch(const MerLengthError& e) {
            std::cerr << e.what() << std::endl;
            continue;
        }
    }
    return retval;
}

bool Filter::join(const Filter& filter) throw(MerLengthError, LowerLevelError) {
    if (filter._lower_level != this->_lower_level) {
        std::ostringstream oss;
        oss << filter._lower_level << " is not " << _lower_level;
        oss << ", Failed joining a filter instance" << std::flush;
        throw LowerLevelError(oss.str());
    }
    if (this->_mer_length == 0) {
        this->_mer_length = filter._mer_length;
    } else if (filter._mer_length != this->_mer_length) {
        std::ostringstream oss;
        oss << filter._mer_length << " is not " << _mer_length;
        oss << ", Failed joining a filter instance";
        throw MerLengthError(oss.str());
    }

    this->_mer_map.insert(filter._mer_map.begin(), filter._mer_map.end());
    return true;
}

std::vector<Filter::score_type> Filter::scores(const Read& read) const {
    std::vector<score_type> retval;
    const int length(read.size() - _mer_length + 1);
    if (_mer_length == 0 || length <= 0) {
        return retval;
    }
    for (int i(0); i < length; i++) {
        const Read sub(read.sub(i, _mer_length));
        if (!sub.isDefinite()) {
            continue;
        }
        retval.push_back(_getScore(sub,_lower_level));
    }
    return retval;
}

bool Filter::check(std::vector<score_type> scores) const {
    if (scores.size() == 0) {
        return false;
    }

    unsigned int lower_total(0), lower_count(0);
    unsigned int upper_total(0), upper_count(0);

    for (std::vector<unsigned int>::const_iterator itr(scores.begin());
            itr != scores.end(); itr++) {
        if (*itr <= _lower_level) {
            lower_count++;
            lower_total += *itr;
        } else {
            upper_count++;
            upper_total += *itr;
            if (lower_count < _lower_interval) {
                lower_count = 0;
            }
        }
    }

    if (upper_count == 0) {
        return true;
    }
    if (lower_count < _lower_interval) {
        return true;
    }

    double upper_average(double(upper_total)/upper_count),
           lower_average(double(lower_total)/lower_count);
    return upper_average < lower_average * _ratio;
}

bool Filter::check(const Read& read) const {
    return check(scores(read));
}

double Filter::average(std::vector<score_type> scores) const {
    if (scores.size() == 0) {
        return 0;
    }

    double total(0.);
    for (std::vector<unsigned int>::const_iterator itr(scores.begin());
            itr != scores.end(); itr++) {
        total += *itr;
    }

    return total/scores.size();
}

double Filter::average(const Read& read) const {
    return average(scores(read));
}

int Filter::_getScore(const Read& read, const score_type default_value) const
        throw(MerLengthError){
    if (read.size() != _mer_length) {
        std::ostringstream oss;
        oss << _mer_length << " is not " << read.size();
        oss << ", Failed getting score of " << read.tostring();
        throw MerLengthError(oss.str());
    }
    int score(0);
    map_type::const_iterator itr(_mer_map.find(read));
    if (itr != _mer_map.end()) {
        score = (*itr).second;
    } else {
        map_type::const_iterator comp(_mer_map.find(read.complement()));
        if (comp != _mer_map.end()) {
            score = (*comp).second;
        } else {
            return default_value;
        }
    }
    return score;
}

} // carl
