#include <vector>
#include <iostream>
#include <string>
#include "../../src/dataio.hpp"
#include "../../src/pairwise_distance_tree.hpp"
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::PairwiseDistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        pstrudel::treeio::write_newick(tree, std::cout);
    }
}
