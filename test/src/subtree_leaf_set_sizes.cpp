#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <unordered_set>
// #include "../../src/platypus-phyloinformary/src/tree.hpp"
#include "../../src/dataio.hpp"
#include "../../src/symmetric_difference_tree.hpp"
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    typedef pstrudel::SymmtericDifferenceTree TreeType;
    auto & out = std::cout;
    std::vector<TreeType> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        tree.calc_subtree_sizes();
        for (auto ndi = tree.postorder_begin(); ndi != tree.postorder_end(); ++ndi) {
            if (ndi.is_leaf()) {
                continue;
            }
            unsigned long nleaves = 0;
            for (auto lfiter = tree.leaf_begin(ndi); lfiter != tree.leaf_end(ndi); ++lfiter) {
                ++nleaves;
            }
            out << ndi->get_label() << "\t" << ndi->get_num_leaves() << "\t" << nleaves << std::endl;
        }
    }
}


