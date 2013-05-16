#include <vector>
#include <iostream>
#include <sstream>
#include <colugo/subprocess.hpp>
#include <platypus/model/treepattern.hpp>
#include <platypus/serialize/newick.hpp>
#include <pstrudel/distancetree.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;
static const unsigned long DEFAULT_NUM_TIPS = 50;
platypus::NewickWriter<pstrudel::DistanceTree> TREE_WRITER;

pstrudel::DistanceTree generate_unbalanced_tree(unsigned long ntips=DEFAULT_NUM_TIPS) {
    pstrudel::DistanceTree tree;
    std::vector<pstrudel::DistanceTree::value_type> leaves;
    leaves.reserve(ntips);
    for (unsigned long i = 0; i < ntips; ++i) {
        leaves.emplace_back();
        auto & leaf = leaves.back();
        leaf.set_label("t" + std::to_string(i+1));
    }
    platypus::build_maximally_unbalanced_tree(tree, leaves.begin(), leaves.end());
    return tree;
}

int check_tree(pstrudel::DistanceTree & tree, int regime) {
    std::string regime_arg;
    if (regime == 0) {
        regime_arg = "mean";
    } else if (regime == 1) {
        regime_arg = "random";
    } else if (regime == 2) {
        regime_arg = "uniform";
    } else if (regime == 3) {
        regime_arg = "anti-coalescent";
    } else {
        throw std::runtime_error("Unsupported coalescent age regime");
    }
    std::ostringstream o;
    TREE_WRITER.write(o, tree);
    colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick"});
    auto result = ps.communicate(o.str());
    // std::cout << ps.returncode() << std::endl;
    if (ps.returncode() != 0) {
    // if (true) {
        std::cerr << result.first;
        std::cerr << result.second;
        return ps.returncode();
    } else {
        return 0;
    }
}

int test_mean_coalescent_ages_on_unbalanced_tree(unsigned long ntips=DEFAULT_NUM_TIPS) {
    auto tree = generate_unbalanced_tree(ntips);
    tree.add_coalescent_edge_lengths(0);
    return check_tree(tree, 0);
}

int main(int, const char * argv[]) {
    TEST_DIR = pstrudel::test::get_test_dir(argv[0]);
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "check-coalescent-edge-lengths.py");
    platypus::bind_standard_interface(TREE_WRITER);
    int fails = 0;
    fails += test_mean_coalescent_ages_on_unbalanced_tree();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
