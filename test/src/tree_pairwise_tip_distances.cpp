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
        tree.calc_pairwise_tip_distances();
        for (auto tip1 = tree.leaf_begin(); tip1 != tree.leaf_end(); ++tip1) {
            for (auto tip2 = tree.leaf_begin(); tip2 != tree.leaf_end(); ++tip2) {
                std::cout << tip1->get_label();
                std::cout << "\t" << tip2->get_label() << "\t";
                std::cout << tree.get_unweighted_pairwise_tip_distance(*tip1, *tip2) << "\t";
                std::cout << std::setprecision(20) << tree.get_weighted_pairwise_tip_distance(*tip1, *tip2) << std::endl;
            }
        }
    }
}


