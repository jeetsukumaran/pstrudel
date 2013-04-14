#include <vector>
#include <iostream>
#include <string>
#include <iomanip>
#include "../../src/dataio.hpp"
#include "../../src/pairwise_distance_tree.hpp"
#include "../../src/split.hpp"
#include "testutils.hpp"

int main(int argc, const char * argv[]) {
    typedef pstrudel::PairwiseDistanceTree TreeType;
    pstrudel::TreeComparisonCalculator<TreeType> tree_calc;
    std::vector<TreeType> trees;
    get_trees(trees, argc, argv);
    unsigned long num_trees = trees.size();
    for (unsigned long tidx1 = 0; tidx1 < num_trees - 1; ++tidx1) {
        auto & tree1 = trees[tidx1];
        for (unsigned long tidx2 = tidx1; tidx2 < num_trees; ++tidx2) {
            auto & tree2 = trees[tidx2];
            auto result = tree_calc.compare_trees(tree1, tree2);
            std::cout << result.false_negatives << "\t";
            std::cout << result.false_positives << "\t";
            std::cout << result.symmetric_difference << "\t";
            std::cout << std::setprecision(20) << result.edge_len_euclidean_distance << "\t";
            std::cout << std::setprecision(20) << result.weighted_rf_distance << "\n";
        }
    }
}


