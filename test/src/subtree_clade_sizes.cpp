#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include <unordered_set>
#include <pstrudel/dataio.hpp>
#include <pstrudel/symmetric_difference_tree.hpp>
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    typedef pstrudel::SymmetricDifferenceTree TreeType;
    auto & out = std::cout;
    std::vector<TreeType> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        tree.calc_subtree_sizes();
        for (auto ndi = tree.postorder_begin(); ndi != tree.postorder_end(); ++ndi) {
            if (ndi.is_leaf()) {
                continue;
            }
            out << ndi->get_label() << "\t" << ndi->get_num_descendent_nodes() << std::endl;
        }
    }
}


