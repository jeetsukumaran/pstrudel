#include <cassert>
#include <vector>
#include <iostream>
#include <sstream>
#include <pstrudel/distancetree.hpp>
#include <platypus/numeric/rng.hpp>
#include "pstrudel_testing.hpp"

int main(int argc, const char * argv[]) {

    int fails = 0;

    typedef std::unordered_multiset<unsigned long> SizesSetType;
    platypus::numeric::RandomNumberGenerator rng;
    std::string test_dir = pstrudel::test::get_test_dir(argv[0]);

    for (unsigned int rep = 0; rep < 20; ++ rep) {
        std::string command = "python " + test_dir + "/calc-multiset-symmetric-difference.py";
        std::string stdout_str = pstrudel::test::execute_external_process(command, true, true);

        std::vector<std::vector<unsigned long>> output;
        std::istringstream s(stdout_str);
        pstrudel::test::read_data_vectors<unsigned long>(s, output );
        assert(output.size() == 3);

        SizesSetType s1;
        s1.insert(output[0].begin(), output[0].end());
        SizesSetType s2;
        s2.insert(output[1].begin(), output[1].end());
        unsigned long expected_diff = output[2][0];

        SizesSetType common;
        SizesSetType unmatched1;
        SizesSetType unmatched2;
        unsigned long obs_diff = pstrudel::DistanceTree::calc_set_symmetric_difference(s1, s2, &common, &unmatched1, &unmatched2);
        fails += pstrudel::test::check_equal(
                expected_diff,
                obs_diff,
                __FILE__,
                __LINE__,
                "test data:\n", stdout_str);
    }

    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
