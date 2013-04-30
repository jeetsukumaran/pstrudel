#include <vector>
#include <iostream>
#include <sstream>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/split.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

int main(int argc, const char * argv[]) {

    typedef pstrudel::DistanceTree  TreeType;

    std::vector<TreeType>  trees;
    std::string test_dir = pstrudel::test::get_test_dir(argv[0]);
    std::string test_data_filepath = pstrudel::test::join_path(test_dir, "data", "trees", "general", "n10-rooted-patterns.nex");
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");

    pstrudel::TreeSplitsDifferences<TreeType> tsd;
    unsigned long num_trees = trees.size();
    unsigned long false_positives;
    unsigned long false_negatives;
    unsigned long symmetric_difference;
    double edge_length_distance;
    std::vector<std::vector<double>> results;
    results.reserve(trees.size() * trees.size());
    for (unsigned int i = 0; i < num_trees - 1; ++i) {
        for (unsigned int j = i + 1; j < num_trees; ++j) {
            results.emplace_back();
            auto & current_result = results.back();
            false_positives = 0;
            false_negatives = 0;
            symmetric_difference = 0;
            edge_length_distance = 0.0;
            tsd.calc_labeled_symmetric_difference(
                    trees[i],
                    trees[j],
                    false_positives,
                    false_negatives,
                    symmetric_difference,
                    edge_length_distance
                    );
            // std::cout << false_positives;
            // std::cout << "\t" << false_negatives;
            // std::cout << "\t" << symmetric_difference;
            // std::cout << "\t" << edge_length_distance;
            // std::cout << std::endl;
            current_result.push_back(false_positives);
            current_result.push_back(false_negatives);
            current_result.push_back(symmetric_difference);
            current_result.push_back(edge_length_distance);
        }
    }

    std::string command = "python " + test_dir + "/calc-splits-distances.py " + test_data_filepath + " nexus";
    std::string reference_results = pstrudel::test::execute_external_process(command, true, true);
    std::istringstream src(reference_results);
    int fails = 0;
    unsigned long row_idx = 0;
    unsigned long col_idx = 0;
    double val = 0.0;
    std::vector<std::string> col_names{"false positives", "false negatives", "symmetric difference", "edge length distance"};
    for (auto & row : results) {
        col_idx = 0;
        for (auto & col : row) {
            src >> val;
            if (!pstrudel::test::is_almost_equal(val, col)) {
                std::cerr << "Comparison " << row_idx + 1 << ": incorrect " << col_names[col_idx] << ": expected " << val << " but found " << col << std::endl;
                fails += 1;
            }
            col_idx += 1;
        }
        row_idx += 1;
    }
    if (fails > 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
