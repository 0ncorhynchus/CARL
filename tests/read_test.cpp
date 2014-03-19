#define BOOST_TEST_MODULE ReadTest

#include <boost/test/included/unit_test.hpp>

#include <string>
#include "../read.hpp"

struct Fixture {
    const std::string sequence_string;
    const std::string complement_string;
    const std::string complement_reverse_string;
    Read read;

    Fixture() :
        sequence_string("tcaggggggttttaatttactttcgtacacagcgtaa"),
        complement_string("agtccccccaaaattaaatgaaagcatgtgtcgcatt"),
        complement_reverse_string("ttacgctgtgtacgaaagtaaattaaaacccccctga"),
        read(sequence_string) {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(constructor) {}

BOOST_AUTO_TEST_CASE(copy_constructor) {
    const Read copy(read);
    BOOST_ASSERT(read == copy);
    BOOST_ASSERT(&read != &copy);
}

BOOST_AUTO_TEST_CASE(copy_assignment) {
    Read copy;
    copy = read;
    BOOST_ASSERT(read == copy);
    BOOST_ASSERT(&read != &copy);
}

BOOST_AUTO_TEST_CASE(tostring) {
    BOOST_CHECK_EQUAL(read.tostring(), sequence_string);
}

BOOST_AUTO_TEST_CASE(size) {
    BOOST_CHECK_EQUAL(read.size(), sequence_string.size());
}

BOOST_AUTO_TEST_CASE(complement) {
    const Read complement(read.complement());
    BOOST_CHECK_EQUAL(complement.size(), read.size());
    BOOST_CHECK_EQUAL(complement.tostring(), complement_string);
}

BOOST_AUTO_TEST_CASE(reverse) {
    const Read complement(read.complement());
    BOOST_CHECK_EQUAL(complement.size(), read.size());
    BOOST_CHECK_EQUAL(complement.tostring(), complement_string);
    const Read complement_reverse(complement.reverse());
    BOOST_CHECK_EQUAL(complement_reverse.size(), complement.size());
    BOOST_CHECK_EQUAL(complement_reverse.tostring(),
            complement_reverse_string);
}

BOOST_AUTO_TEST_CASE(getBaseAt) {
    char bases[] = {'a', 'c', 'g', 't'};
    for (int i = 0; i < read.size(); i++) {
        const int at(read.getBaseAt(i));
        BOOST_ASSERT(at >= 0 && at < 4);
        BOOST_CHECK_EQUAL(bases[at], sequence_string.at(i));
    }
}

BOOST_AUTO_TEST_CASE(isDefinite) {
    BOOST_CHECK(read.isDefinite());
    const Read invalid_sequence("atgcugatc");
    BOOST_CHECK(!invalid_sequence.isDefinite());
}

BOOST_AUTO_TEST_SUITE_END()
