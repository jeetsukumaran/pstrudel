#include <vector>
#include <iostream>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>
#include "supportlib.hpp"

int main() {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    pstrudel::treeio::read_from_filepath(trees, "data/trees/basic-patterns/n06-rooted-patterns.nex", "nexus");
    if (trees.size() != 6) {
        return pstrudel::test::fail_test(__FILE__,
                6, trees.size(), "Incorrect number of trees read from file");
    }
    for (auto & tree : trees) {
        if (tree.get_num_tips() != 6) {
            return pstrudel::test::fail_test(__FILE__,
                    6, tree.get_num_tips(), "Incorrect number of tips recorded");
        }
    }
    return EXIT_SUCCESS;
}
