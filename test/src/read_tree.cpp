#include <vector>
#include <iostream>
#include <string>
#include "testutils.hpp"
#include <pstrudel/dataio.hpp>
#include <pstrudel/pairwise_distance_tree.hpp>

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::PairwiseDistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        pstrudel::treeio::write_newick(tree, std::cout);
    }
}
