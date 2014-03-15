#define BOOST_TEST_MODULE FastaTest

#include <boost/test/included/unit_test.hpp>

#include <fstream>
#include <sstream>
#include <algorithm>
#include "../fasta.hpp"

struct Fixture {
	const std::string filename;
	Fasta fasta;

	Fixture() :
		filename("sample.fasta"),
		fasta(filename)
	{
	}

};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(constructor) {}

BOOST_AUTO_TEST_CASE(item) {
	const std::string info("information");
	const std::string sequence("agtcgtacgtcgggtcccaaagtgagtgt");
	Fasta::Item item(info, sequence);
	BOOST_CHECK_EQUAL(item.getInfo(), info);
	BOOST_CHECK_EQUAL(item.getRead().tostring(), sequence);
}

BOOST_AUTO_TEST_CASE(getItemStrings) {
	std::pair<std::string, std::string> item_string(
			fasta.getItemStrings());
	std::ifstream ifs(filename);
	std::string buff;
	getline(ifs, buff);
	std::ostringstream oss;
	oss << ">" << item_string.first;
	BOOST_CHECK_EQUAL(oss.str(), buff);
	getline(ifs, buff);
	BOOST_CHECK_EQUAL(item_string.second, buff);
	ifs.close();
}

BOOST_AUTO_TEST_CASE(getItem) {
	Fasta::Item item(fasta.getItem());
	std::ifstream ifs(filename);
	std::string buff;
	getline(ifs, buff);
	std::ostringstream oss;
	oss << ">" << item.getInfo();
	BOOST_CHECK_EQUAL(oss.str(), buff);
	getline(ifs, buff);
	transform(buff.begin(), buff.end(), buff.begin(), tolower);
	BOOST_CHECK_EQUAL(item.getRead().tostring(), buff);
	ifs.close();
}

BOOST_AUTO_TEST_SUITE_END()
