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

        // Expected script output is 3 lines
        // 1: vector of numbers
        // 2: vector of numbers
        // 3: difference between 1 and 2, i.e. number of elements remaining
        // once common elements are removed based on matches.
        // E.g.
        //
        //      74 15 51 38 65 1 68 57 69 54 25 46 9 64 48 25 29 85 33 96 94 42 18 37 95 31 43 100 16 59 3 32 26 59 10 19 74 3 32 27 53 93 30 76 75 80 86 71 11 6 35 65 63 75 91 93 74 66 55 46 21 79 12 28 9 26 36 85 21 85 55 56 19 96 14 24 50 3 97 5 32 36 56 90 52 89 13 38 85 63 8 58 23 27 29 47 43 41 2 60 67 73 17 100 39 22 68 80 16 52 1 81 81 38 29 60 49 11 64
        //      73 96 48 64 77 52 91 96 86 64 39 50 93 80 83 53 22 44 44 47 70 44 86 51 98 63 34 86 78 72 83 30 97 95 66 20 17 94 86 78 9 14 8 98 31 72 75 33 49 40 50 52 58 60 60 7 60 59 78 22 61 40 37 71 44 61 77 9 82 85 56 2 38 2 61 65 58 40 66 54 31 96 33 77 39 20 94 95 82 87 99 18 60 10 15 52 14 57 41 94
        //      119
        //
        // E.g.
        //      1 1 2 2 3 3 6 7 8
        //      1 2 2 3 4 5 6 8 9
        //      common = 1, 2, 2, 3, 6, 8
        //      remaining1 = 3, 6, 7
        //      remaining2 = 3, 4, 9
        //      diff = 6
        //
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
