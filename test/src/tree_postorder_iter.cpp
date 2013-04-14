#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include "../../src/dataio.hpp"
#include "../../src/pairwise_distance_tree.hpp"
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    std::vector<pstrudel::PairwiseDistanceTree> trees;
    get_trees(trees, argc, argv);
    for (auto & tree : trees) {
        int postorder_count = 0;
        for (auto ndi = tree.postorder_begin(); ndi != tree.postorder_end(); ++ndi, ++postorder_count) {
            std::cout << ndi->get_label() << "\t" << std::setprecision(8) << ndi->get_edge_length() << std::endl;
        }
    }
}


