#define BOOST_TEST_MODULE FilterTest

#include <boost/test/included/unit_test.hpp>

#include "../filter.hpp"

struct Fixture {
    const std::string filename, countname;
    Filter filter;

    Fixture() :
        filename("samples/sample.fasta"),
        countname("samples/sample.count"),
        filter()
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(constructor) {}

bool check_exception(const Filter::MerLengthError& e) {
    std::cout << e.what() << std::endl;
    return true;
}

BOOST_AUTO_TEST_CASE(insertMer) {
    Read read("acgttttgggaacgcgcgttgtgtaact");
    Read invalid_read("acgtgnactggtgggc");
    int score(52);
    try {
        BOOST_CHECK(!filter.insertMer(invalid_read, score));
        BOOST_CHECK(filter.insertMer(read, score));
    } catch(const Filter::MerLengthError& e) {
        BOOST_FAIL(e.what());
    }
    Read length_error("acgtggggggttttaccactgt");
    BOOST_CHECK_EXCEPTION(filter.insertMer(length_error, score),
            Filter::MerLengthError, check_exception);
}

BOOST_AUTO_TEST_CASE(insertMers) {
    filter = Filter(10,20,2.);
    Fasta count(countname);
    BOOST_CHECK(filter.insertMers(count));
}

BOOST_AUTO_TEST_CASE(size) {
    BOOST_CHECK_EQUAL(filter.size(), 0);
    Read read("acgtgactgta");
    int score(52);
    try {
        BOOST_CHECK(filter.insertMer(read, score));
    } catch(const Filter::MerLengthError& e) {
        BOOST_FAIL(e.what());
    }
    BOOST_CHECK_EQUAL(filter.size(), 1);
}

BOOST_AUTO_TEST_CASE(join) {
    Filter other;
    Fasta fasta(countname);
    BOOST_CHECK(other.insertMers(fasta));
    try {
        BOOST_CHECK(filter.join(other));
    } catch(const Filter::MerLengthError& e) {
        BOOST_FAIL(e.what());
    } catch(const Filter::LowerLevelError& e) {
        BOOST_FAIL(e.what());
    }
    BOOST_CHECK_EQUAL(filter.size(), other.size());
}

BOOST_AUTO_TEST_CASE(scores) {
    filter = Filter(10,20,2.);
    Fasta count(countname);
    BOOST_CHECK(filter.insertMers(count));
    Fasta fasta(filename);
    Fasta::Item item(fasta.getItem());
    std::vector<unsigned int> scores(filter.scores(item.getRead()));
    bool flg(false);
    for (std::vector<unsigned int>::const_iterator itr(scores.begin());
            itr != scores.end(); itr++) {
        std::cout << *itr << " ";
        if (*itr > 10)
            flg = true;
    }
    std::cout << std::endl;
    BOOST_CHECK(flg);
}

BOOST_AUTO_TEST_CASE(check) {
    filter = Filter(10,20,2.);
    std::vector<unsigned int> scores;
    for (int i(0); i < 100; i++) {
        scores.push_back(11);
    }
    BOOST_CHECK(filter.check(scores));
    scores.clear();
    for (int i(0); i < 100; i++) {
        scores.push_back(10);
    }
    BOOST_CHECK(filter.check(scores));
    scores.clear();
    for (int i(0); i < 100; i++) {
        if (i >= 10 && i < 31)
            scores.push_back(10);
        else
            scores.push_back(40);
    }
    BOOST_CHECK(!filter.check(scores));
    scores.clear();
    for (int i(0); i < 100; i++) {
        if (i >= 10 && i < 31)
            scores.push_back(10);
        else
            scores.push_back(40);
    }
    BOOST_CHECK(!filter.check(scores));
}

BOOST_AUTO_TEST_CASE(average) {
    filter = Filter(1,20,2.);
    std::vector<unsigned int> scores;
    for (int i(0); i < 100; i++) {
        scores.push_back(i);
    }
    BOOST_CHECK_EQUAL(filter.average(scores), 49.5);
}

BOOST_AUTO_TEST_SUITE_END()
