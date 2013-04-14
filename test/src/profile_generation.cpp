#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "../../src/profile.hpp"

int main(int argc, const char * argv[]) {
    std::vector<double> data;
    data.insert(data.end(), std::istream_iterator<double>(std::cin), std::istream_iterator<double>());
    unsigned long num_points = 1000;
    if (argc > 1) {
        num_points = std::stol(argv[1]);
    }
    pstrudel::Profile profile;
    profile.add_subprofile_source("test", data.begin(), data.end(), num_points);
    profile.build_subprofiles();
    auto profile_data = profile.get_subprofile("test");
    std::cout << std::setprecision(20);
    std::copy(profile_data.cbegin(), profile_data.cend(), std::ostream_iterator<double>(std::cout, "\n"));
}
