#include <vector>
#include <iostream>
#include <pstrudel/profile.hpp>
#include "pstrudel_testing.hpp"

int test_profile_distance0() {
    int fails=0;
    pstrudel::Profile   profile1(0);
    pstrudel::Profile   profile2(0);
    std::vector<double> raw_data{3, 7, 9, 10, 19, 17, 19, 1};
    profile1.set_data(raw_data.begin(), raw_data.end());
    profile2.set_data(raw_data.begin(), raw_data.end());
    fails += platypus::testing::compare_equal(
            0,
            profile1.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_equal(
            0,
            profile2.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_equal(
            0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    profile1.fix_size(1000);
    fails += platypus::testing::compare_equal(
            0,
            profile1.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_equal(
            0,
            profile2.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_equal(
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
    fails += platypus::testing::compare_equal(
            0.0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_equal(
            0.0,
            profile2.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    profile1.fix_size(1000);
    // std::cout << "***" << std::endl;
    // std::cout << "***" << std::endl;
    // pstrudel::test::write_container(std::cout, profile1.get_profile(1000));
    // std::cout << std::endl;
    // std::cout << "***" << std::endl;
    // std::cout << std::endl;
    // pstrudel::test::write_container(std::cout, profile2.get_profile(1000));
    // std::cout << std::endl;
    // std::cout << "***" << std::endl;
    // std::cout << "***" << std::endl;
    fails += platypus::testing::compare_almost_equal(
            0.0,
            profile1.get_distance(profile2),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    fails += platypus::testing::compare_almost_equal(
            0.0,
            profile2.get_distance(profile1),
            __FILE__,
            __LINE__,
            "incorrect profile distance");
    return fails;
}

// options(digits=22)
// euc.dist = function(x1, x2) { sqrt(sum((x1 - x2) ^ 2))}
// a = c(0, 1.5, 3, 5, 7, 8, 9, 9.5, 10)
// b = c(0, 3, 6, 7.5, 9, 11, 13, 15, 17)
// euc.dist(a, b)
// 11.2138307460029018614
int test_profile_distance2() {

    int fails=0;

    std::vector<double> raw_data1{3, 7, 9, 10};
    pstrudel::Profile   profile1(0);
    profile1.set_data(raw_data1.begin(), raw_data1.end());

    std::vector<double> raw_data2{6, 9, 13, 17};
    pstrudel::Profile   profile2(0);
    profile2.set_data(raw_data2.begin(), raw_data2.end());

    fails += platypus::testing::compare_almost_equal(
            8.831760866327847736557,
            profile1.get_distance(profile2, false),
            __FILE__,
            __LINE__,
            "incorrect profile distance");

    std::vector<double> raw_data3{1.5, 3, 5, 7, 8, 9, 9.5, 10};
    pstrudel::Profile   profile3(0);
    profile3.set_data(raw_data3.begin(), raw_data3.end());
    // std::cout << "***" << std::endl;
    // pstrudel::test::write_container(std::cout, profile3.get_profile());
    // std::cout << std::endl;
    // pstrudel::test::write_container(std::cout, profile2.get_profile(profile3.data_size()));
    // std::cout << std::endl;
    // std::cout << "***" << std::endl;
    fails += platypus::testing::compare_almost_equal(
            0.0,
            profile1.get_distance(profile3, false),
            __FILE__,
            __LINE__,
            "incorrect profile distance");

    fails += platypus::testing::compare_almost_equal(
            11.2138307460029018614,
            profile2.get_distance(profile3, false),
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
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
