#include <vector>
#include <iostream>
#include <pstrudel/treeshape.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

int main() {
    typedef pstrudel::TreeShape  TreeType;
    std::vector<TreeType>  trees;
    pstrudel::treeio::read_from_filepath(trees, "data/trees/general/n06-rooted-patterns.nex", "nexus");
    int fails = 0;
    if (trees.size() != 6) {
        fails += platypus::testing::compare_equal(
                6UL,
                trees.size(),
                __FILE__,
                __LINE__,
                "Incorrect number of trees read from file");
    }
    for (auto & tree : trees) {
        if (tree.get_num_tips() != 6) {
            fails += platypus::testing::compare_equal(
                    6UL,
                    tree.get_num_tips(),
                    __FILE__,
                    __LINE__,
                    "Incorrect number of tips recorded");
        }
    }
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
