#include <vector>
#include <string>
#include <iostream>
#include <platypus/parse/newick.hpp>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

int main(int argc, const char * argv[]) {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    std::string test_dir = pstrudel::test::get_test_dir(argv[0]);
    std::string test_data_filepath = pstrudel::test::join_path(test_dir, "data", "trees", "general", "n10-rooted-patterns.nex");
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");

    std::map<unsigned long, std::map<unsigned long, double>> uwt_dists;
    std::map<unsigned long, std::map<unsigned long, double>> wt_dists;
    std::set<std::pair<unsigned long, unsigned long>> comparisons;
    for (unsigned tidx1 = 0; tidx1 < trees.size(); ++tidx1) {
        for (unsigned tidx2 = 0; tidx2 < trees.size(); ++tidx2) {
            auto & tree1 = trees[tidx1];
            auto & tree2 = trees[tidx2];
            uwt_dists[tidx1][tidx2] = tree1.get_unweighted_pairwise_tip_profile_distance(tree2);
            wt_dists[tidx1][tidx2] = tree1.get_weighted_pairwise_tip_profile_distance(tree2);
            comparisons.insert(std::make_pair(tidx1, tidx2));
        }
    }

    std::string command = "python " + test_dir + "/calc-tree-pairwise-profile-distance.py " + test_data_filepath + " nexus";
    std::string reference_results = pstrudel::test::execute_external_process(command, true, true);
    std::istringstream src(reference_results);
    int fails = 0;
    unsigned long tidx1 = 0;
    unsigned long tidx2 = 0;
    double exp_uwt_d = 0.0;
    double exp_wt_d = 0.0;
    while (src.good()) {
        src >> tidx1;
        src >> tidx2;
        src >> exp_uwt_d;
        src >> exp_wt_d;
        // std::cerr << tidx1 << "\t" << tidx2 << "\t" << exp_uwt_d << "\t" << exp_wt_d << std::endl;
        fails += platypus::testing::compare_equal(
                exp_uwt_d,
                uwt_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect unweighted pairwise tip profile distance for tree ", tidx1, " and tree ", tidx2);
        fails += platypus::testing::compare_equal(
                exp_wt_d,
                wt_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect weighted pairwise tip profile distance for tree ", tidx1, " and tree ", tidx2);
        comparisons.erase(std::make_pair(tidx1, tidx2));
    }
    fails += platypus::testing::compare_equal(
            0UL,
            comparisons.size(),
            __FILE__,
            __LINE__,
            "Not all trees checked");

    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }

}
