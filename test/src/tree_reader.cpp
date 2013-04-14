#include <vector>
#include <iostream>
#include <string>
#include "../../src/platypus-phyloinformary/src/ncl_reader.hpp"
#include "../../src/dataio.hpp"
#include "../../src/pairwise_distance_tree.hpp"
#include "testutils.hpp"

const char * DATA_FORMAT = "nexus";

template <typename T>
void test_raw_tree_reader(std::istream& src) {
    std::vector<T> trees;
    std::function<T& ()> tree_factory ( [&trees] () -> T& { trees.emplace_back(); return trees.back(); } );
    platypus::NclTreeReader<T> tree_reader(tree_factory);
    tree_reader.set_node_value_label_func(&T::value_type::set_label);
    tree_reader.set_node_value_edge_length_func(&T::value_type::set_edge_length);
    tree_reader.read_from_stream(src, DATA_FORMAT);
    for (auto & tree : trees) {
        pstrudel::treeio::write_newick(tree, std::cout);
    }
}

int main(int argc, const char * argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <TREEFILE> <CODE-BRANCH>" << std::endl;
        exit(1);
    }
    std::ifstream src(argv[1]);
    unsigned long code_branch = std::stol(argv[2]);;

    // Tests:
    // 0    -   TreeReader, nodes: std::allocator
    // 1    -   TreeVectorReader, nodes: std::allocator
    // 2    -   TreePointerVectorReader, nodes: std::allocator
    // 3    -   TreePointerVectorReader, nodes: custom

    if (code_branch == 0) {
        test_raw_tree_reader<pstrudel::PairwiseDistanceTree>(src);
    // } else if (code_branch == 1) {
    // } else if (code_branch == 2) {
    // } else if (code_branch == 3) {
    // } else if (code_branch == 4) {
    // } else if (code_branch == 5) {
    // } else if (code_branch == 6) {
    // } else if (code_branch == 7) {
    // } else if (code_branch == 8) {
    // } else if (code_branch == 9) {
    } else {
        std::cerr << "Invalid code branch: " << code_branch << std::endl;
        exit(1);
    }
}
