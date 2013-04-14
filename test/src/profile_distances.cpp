#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "testutils.hpp"
#include "../../src/colugo-utilities/src/utility.hpp"
#include "../../src/colugo-utilities/src/textutil.hpp"
#include "../../src/profile.hpp"

int main(int argc, const char * argv[]) {
    unsigned long num_subprofiles = 1;
    if (argc > 1) {
        num_subprofiles = std::stol(argv[1]);
    }

    std::vector<std::vector<double>> data;
    read_data_vectors(std::cin, data);

    unsigned long expected_num_profiles = static_cast<unsigned long>(data.size() / num_subprofiles);
    unsigned long expected_num_subprofiles = static_cast<unsigned long>(data.size() / expected_num_profiles);
    unsigned long expected_num_rows = expected_num_profiles * expected_num_subprofiles;
    if (expected_num_rows != data.size()) {
        colugo::colugo_abort("Expecting ",
                expected_num_rows,
                " of data rows for ",
                expected_num_profiles,
                " profiles with ",
                num_subprofiles,
                " subprofiles per profile, but instead found ",
                data.size());
    }

    std::vector<pstrudel::Profile> profiles;
    unsigned long data_row_idx = 0;
    for (unsigned long profile_idx = 0; profile_idx < expected_num_profiles; ++profile_idx) {
        profiles.emplace_back();
        auto & profile = profiles.back();
        for (unsigned long subprofile_idx = 0; subprofile_idx < num_subprofiles; ++subprofile_idx) {
            std::string label = colugo::textutil::join("-", "subprofile", subprofile_idx+1);
            // std::cerr << "(" << data_row_idx+1 << " of " << data.size() << ") ";
            // std::cerr << profile_idx << ":" << subprofile_idx << ":  '" << label << "'" << std::endl;
            profile.insert_subprofile(label.c_str(), data[data_row_idx].cbegin(), data[data_row_idx].cend());
            ++data_row_idx;
        }
    }

    for (unsigned long idx1 = 0; idx1 < profiles.size() - 1; ++idx1) {
        for (unsigned long idx2 = idx1 + 1; idx2 < profiles.size(); ++idx2) {
            std::cout << std::setprecision(20) << profiles[idx1].calc_distance(profiles[idx2]) << std::endl;
        }
    }
}
