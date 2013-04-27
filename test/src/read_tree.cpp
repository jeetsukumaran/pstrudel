#include <vector>
#include <iostream>
#include <string>
#include "testutils.hpp"
#include <pstrudel/dataio.hpp>
#include <pstrudel/distancetree.hpp>

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::DistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        pstrudel::treeio::write_newick(tree, std::cout);
    }
}
