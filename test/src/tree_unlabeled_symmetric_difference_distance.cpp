#include <vector>
#include <string>
#include <iostream>
#include <platypus/parse/newick.hpp>
#include <pstrudel/treeshape.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

int main(int argc, const char * argv[]) {
    typedef pstrudel::TreeShape  TreeType;
    std::vector<TreeType>  trees;
    std::string test_dir = pstrudel::test::get_test_dir(argv[0]);
    std::string test_data_filepath = pstrudel::test::join_path(test_dir, "data", "trees", "general", "n10-rooted-patterns.nex");
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");

    std::map<unsigned long, std::map<unsigned long, double>> obs_dists;
    std::set<std::pair<unsigned long, unsigned long>> comparisons;
    for (unsigned tidx1 = 0; tidx1 < trees.size(); ++tidx1) {
        for (unsigned tidx2 = 0; tidx2 < trees.size(); ++tidx2) {
            auto & tree1 = trees[tidx1];
            auto & tree2 = trees[tidx2];
            obs_dists[tidx1][tidx2] = tree1.get_unlabeled_symmetric_difference(tree2);
            comparisons.insert(std::make_pair(tidx1, tidx2));
        }
    }
    std::string command = "python " + test_dir + "/calc-tree-unlabeled-symmetric-difference.py " + test_data_filepath + " nexus";
    std::string reference_results = pstrudel::test::execute_external_process(command, true, true);
    std::istringstream src(reference_results);
    int fails = 0;
    unsigned long tidx1 = 0;
    unsigned long tidx2 = 0;
    unsigned long exp_dist = 0;
    while (src.good()) {
        src >> tidx1;
        src >> tidx2;
        src >> exp_dist;
        fails += platypus::testing::compare_equal(
                exp_dist,
                obs_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect unlabeled symmetric difference distance for tree ", tidx1, " and tree ", tidx2);
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
