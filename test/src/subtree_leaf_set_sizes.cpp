#include <set>
#include <set>
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

    std::map<unsigned long, std::map<std::string, unsigned long>> obs_subtree_leaf_set_sizes;
    std::map<unsigned long, std::multiset<std::string>> obs_subtree_nodes;
    unsigned long tree_idx = 0;
    for (auto & tree : trees) {
        tree.get_symmetric_difference_calculator().calc_subtree_leaf_set_sizes();
        for (auto nd_iter = tree.postorder_begin(); nd_iter != tree.postorder_end(); ++nd_iter) {
            if (nd_iter.is_leaf()) {
                continue;
            }
            auto & nd = *nd_iter;
            obs_subtree_leaf_set_sizes[tree_idx][nd.get_label()] = nd.get_num_leaves();
            obs_subtree_nodes[tree_idx].insert(nd.get_label());
        }
        tree_idx += 1;
    }

    std::string command = "python " + test_dir + "/calc-subtree-leaf-set-sizes.py " + test_data_filepath + " nexus";
    std::string reference_results = pstrudel::test::execute_external_process(command, true, true);
    std::istringstream src(reference_results);
    std::map<unsigned long, std::map<std::string, unsigned long>> exp_subtree_leaf_set_sizes;
    unsigned long script_tree_idx;
    std::string script_nd_label;
    unsigned long script_leaf_set_size;
    while (src.good()) {
        src >> script_tree_idx;
        src >> script_nd_label;
        src >> script_leaf_set_size;
        exp_subtree_leaf_set_sizes[script_tree_idx][script_nd_label] = script_leaf_set_size;
    }

    int fails = 0;

    fails += platypus::testing::compare_equal( exp_subtree_leaf_set_sizes.size(), obs_subtree_leaf_set_sizes.size(), __FILE__, __LINE__, "Mismatch in number of trees");

    for (auto & exp_tree_iter : exp_subtree_leaf_set_sizes) {
        auto & tree_idx = exp_tree_iter.first;
        auto & exp_lss = exp_tree_iter.second;
        if (obs_subtree_leaf_set_sizes.find(tree_idx) == obs_subtree_leaf_set_sizes.end()) {
            fails += platypus::testing::fail_test(__FILE__, __LINE__, tree_idx, "(not found)", "Tree index not found");
        } else {
            auto obs_lss = obs_subtree_leaf_set_sizes[tree_idx];
            for (auto & exp_lss_iter : exp_lss) {
                auto & nd_label = exp_lss_iter.first;
                auto & exp_size = exp_lss_iter.second;
                auto obs_lss_iter = obs_lss.find(nd_label);
                if (obs_lss_iter == obs_lss.end()) {
                    fails += platypus::testing::fail_test(__FILE__, __LINE__, nd_label, "(not found)", "Tree ", tree_idx, ": '", nd_label, "' not found");
                } else {
                    auto obs_size = obs_lss_iter->second;
                    fails += platypus::testing::compare_equal(
                            exp_size,
                            obs_size,
                            __FILE__,
                            __LINE__,
                            "Tree ", tree_idx, ": '", nd_label, "' subtree leaf set size");
                    obs_subtree_leaf_set_sizes[tree_idx].erase(obs_subtree_leaf_set_sizes[tree_idx].find(nd_label));
                }
            }
            fails += platypus::testing::compare_equal(
                    0UL,
                    obs_subtree_leaf_set_sizes[tree_idx].size(),
                    __FILE__,
                    __LINE__,
                    "Tree ", tree_idx, ": nodes remaining");
        }
    }
    if (fails == 0) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}
