#include <vector>
#include <iostream>
#include <string>
#include <pstrudel/distancetree.hpp>
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::DistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & original_tree : trees) {
        pstrudel::DistanceTree clone(original_tree);
    }
}
