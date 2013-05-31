#include <map>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <pstrudel/treeshape.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

int main(int argc, const char * argv[]) {

    typedef pstrudel::TreeShape  TreeType;

    std::vector<TreeType>  trees;
    std::string test_dir = pstrudel::test::get_test_dir(argv[0]);
    std::string test_data_filepath = pstrudel::test::join_path(test_dir, "data", "trees", "general", "n10-rooted-patterns.nex");
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");

    for (auto & tree : trees) {
        std::map< std::pair<std::string, std::string>, unsigned long >  uw_pwdists;
        std::map< std::pair<std::string, std::string>, double >         we_pwdists;
        unsigned long anon_nd_idx = 0;
        std::string anon_nd_prefix = "z";
        auto f = [&uw_pwdists, &we_pwdists, &anon_nd_idx, &anon_nd_prefix] (
                pstrudel::DistanceNodeValue & n1,
                pstrudel::DistanceNodeValue & n2,
                unsigned long u,
                double w) {
            if (n1.get_label().empty()) {
                n1.set_label(anon_nd_prefix + std::to_string(++anon_nd_idx));
            }
            if (n2.get_label().empty()) {
                n2.set_label(anon_nd_prefix + std::to_string(++anon_nd_idx));
            }
            auto key1 = std::make_pair(n1.get_label(), n2.get_label());
            auto key2 = std::make_pair(n2.get_label(), n1.get_label());
            uw_pwdists[key1] = u;
            uw_pwdists[key2] = u;
            we_pwdists[key1] = w;
            we_pwdists[key2] = w;
        };
        tree.calc_pairwise_tip_distances(f);
        std::vector<std::string> leaves;
        for (auto ndi = tree.leaf_begin(); ndi != tree.leaf_end(); ++ndi) {
            leaves.push_back(ndi->get_label());
        }
        std::string command = "python " + test_dir + "/calc-pairwise-tip-distances.py " + test_data_filepath + " nexus";
        std::string reference_results = pstrudel::test::execute_external_process(command, true, true);
        std::istringstream src(reference_results);
        int fails = 0;
        std::string tip1;
        std::string tip2;
        unsigned long r_uwdist;
        double r_wedist;
        std::sort(leaves.begin(), leaves.end());
        for (unsigned long idx1 = 0; idx1 < leaves.size() - 1; ++idx1) {
            for (unsigned long idx2= idx1 + 1; idx2 < leaves.size(); ++idx2) {
                if (!src.good()) {
                    throw std::runtime_error("Unexpected end of source");
                }
                src >> tip1;
                src >> tip2;
                src >> r_uwdist;
                src >> r_wedist;
                auto key = std::make_pair(tip1, tip2);
                fails += platypus::testing::compare_equal(
                        r_uwdist,
                        uw_pwdists[key],
                        __FILE__,
                        __LINE__,
                        "Incorrect unweighted pairwise tip distance for '", tip1, "' and '", tip2, "'");
                fails += platypus::testing::compare_equal(
                        r_wedist,
                        we_pwdists[key],
                        __FILE__,
                        __LINE__,
                        "Incorrect unweighted pairwise tip distance for '", tip1, "' and '", tip2, "'");
            }
        }
        if (fails != 0) {
            return EXIT_FAILURE;
        } else {
            return EXIT_SUCCESS;
        }
    }
}
