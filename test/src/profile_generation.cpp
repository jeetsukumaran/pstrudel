#include <vector>
#include <iostream>
#include <pstrudel/profile.hpp>
#include "pstrudel_testing.hpp"

int test_basic_profile_population() {
    int fails=0;
    pstrudel::Profile   profile(0);
    std::vector<double> raw_data{3, 7, 9, 10};
    profile.set_data(raw_data.begin(), raw_data.end());
    raw_data.insert(raw_data.begin(), 0); // profile generation adds '0'
    fails += pstrudel::test::check_equal(
            raw_data,
            profile.get_profile(raw_data.size()),
            __FILE__,
            __LINE__,
            "incorrect raw-data sized profile");
    return fails;
}

int test_basic_profile_interpolation() {
    int fails=0;
    pstrudel::Profile   profile(0);
    std::vector<double> raw_data{3, 7, 9, 10};
    profile.set_data(raw_data.begin(), raw_data.end());
    std::vector<double> expected{0, 1.5, 3, 5, 7, 8, 9, 9.5, 10};
    fails += pstrudel::test::check_equal(
            expected,
            profile.get_profile(raw_data.size() * 2),
            __FILE__,
            __LINE__,
            "incorrect interpolated profile");
    return fails;
}

int main() {
    int fails = 0;
    fails += test_basic_profile_population();
    fails += test_basic_profile_interpolation();
    if (fails > 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
