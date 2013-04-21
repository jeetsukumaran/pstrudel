#include <vector>
#include <iostream>
#include <string>
#include <pstrudel/dataio.hpp>
#include <pstrudel/pairwise_distance_tree.hpp>
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::PairwiseDistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        int leaf_count = 0;
        for (auto ndi = tree.leaf_begin(); ndi != tree.leaf_end(); ++ndi, ++leaf_count) {
            std::cout << leaf_count + 1 << "\t" << ndi->get_label() << std::endl;
            if (leaf_count > 40) {
                std::cerr << "terminating: too many nodes" << std::endl;
                exit(1);
            }
        }
    }
}
