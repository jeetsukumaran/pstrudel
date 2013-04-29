#include <vector>
#include <iostream>
#include <pstrudel/profile.hpp>
#include "supportlib.hpp"

int test_profile_distance0() {
    int fails=0;
    pstrudel::Profile   profile1(0);
    pstrudel::Profile   profile2(0);
    std::vector<double> raw_data{3, 7, 9, 10, 19, 17, 19, 1};
    profile1.set_data(raw_data.begin(), raw_data.end());
    profile2.set_data(raw_data.begin(), raw_data.end());
    fails += pstrudel::test::check_equal(
            0,
            profile1.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_equal(
            0,
            profile2.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_equal(
            0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    profile1.fix_size(1000);
    fails += pstrudel::test::check_equal(
            0,
            profile1.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_equal(
            0,
            profile2.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_equal(
            0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    return fails;
}

int test_profile_distance1() {
    int fails=0;
    pstrudel::Profile   profile1(0);
    pstrudel::Profile   profile2(0);
    std::vector<double> raw_data1{3, 7, 9, 10};
    profile1.set_data(raw_data1.begin(), raw_data1.end());
    std::vector<double> raw_data2{1.5, 3, 5, 7, 8, 9, 9.5, 10};
    profile2.set_data(raw_data2.begin(), raw_data2.end());
    fails += pstrudel::test::check_equal(
            0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_equal(
            0,
            profile2.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    profile1.fix_size(1000);
    fails += pstrudel::test::check_almost_equal(
            0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += pstrudel::test::check_almost_equal(
            0,
            profile2.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    return fails;
}

int test_profile_distance2() {
    int fails=0;
    pstrudel::Profile   profile1(0);
    pstrudel::Profile   profile2(0);
    std::vector<double> raw_data1{3, 7, 9, 10};
    profile1.set_data(raw_data1.begin(), raw_data1.end());
    std::vector<double> raw_data2{6, 9, 13, 17};
    profile2.set_data(raw_data2.begin(), raw_data2.end());
    double dist = 8.83176;
    fails += pstrudel::test::check_almost_equal(
            dist,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    return fails;
}


int main() {
    int fails = 0;
    fails += test_profile_distance0();
    fails += test_profile_distance1();
    fails += test_profile_distance2();
    if (fails > 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
