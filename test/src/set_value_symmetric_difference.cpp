#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <colugo/utility.hpp>
#include <colugo/textutil.hpp>
#include <pstrudel/distancetree.hpp>
#include "testutils.hpp"

void dump_results(pstrudel::DistanceTree::SizesSetType & s) {
    std::vector<unsigned long> v(s.cbegin(), s.cend());
    std::sort(v.begin(), v.end());
    std::copy(v.cbegin(), v.cend(), std::ostream_iterator<unsigned long>(std::cout, " "));
}

int main(int argc, const char * argv[]) {
    std::vector<std::vector<unsigned long>> data;
    read_data_vectors(std::cin, data);
    if (data.size() != 2) {
        colugo::colugo_abort("Expecting exactly 2 rows of data, but found ", data.size());
    }
    pstrudel::DistanceTree::SizesSetType set1(data[0].begin(), data[0].end());
    pstrudel::DistanceTree::SizesSetType set2(data[1].begin(), data[1].end());
    pstrudel::DistanceTree::SizesSetType common;
    pstrudel::DistanceTree::SizesSetType unmatched1;
    pstrudel::DistanceTree::SizesSetType unmatched2;
    unsigned long sdiff = pstrudel::DistanceTree::calc_set_symmetric_difference(set1, set2, &common, &unmatched1, &unmatched2);
    std::cout << sdiff << std::endl;

    dump_results(common);
    std::cout << std::endl;

    pstrudel::DistanceTree::SizesSetType combined_unmatched;
    combined_unmatched.insert(unmatched1.cbegin(), unmatched1.cend());
    combined_unmatched.insert(unmatched2.cbegin(), unmatched2.cend());
    dump_results(combined_unmatched);
    std::cout << std::endl;
    std::cout << combined_unmatched.size() << std::endl;
}
