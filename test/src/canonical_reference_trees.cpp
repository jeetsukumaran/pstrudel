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

int check_tree(pstrudel::DistanceTree & tree, int regime, const std::string & label) {
    std::string regime_arg;
    if (regime == 0) {
        regime_arg = "mean-coalescent";
    } else if (regime == 1) {
        regime_arg = "converse-coalescent";
    } else if (regime == 2) {
        regime_arg = "uniform";
    } else {
        throw std::runtime_error("Unsupported coalescent age regime");
    }
    std::ostringstream o;
    TREE_WRITER.write(o, tree);
    colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick", "-l", label, "-r", regime_arg});
    auto result = ps.communicate(o.str());
    if (ps.returncode() != 0) {
        std::cerr << "(test '" << label << "' returned error: " << ps.returncode() << ")\n";
        TREE_WRITER.write(std::cerr, tree);
        std::cerr << std::endl;
        std::cerr << result.first;
        std::cerr << result.second;
        return 1;
    } else {
        return 0;
    }
}

int test_mean_coalescent_ages_on_unbalanced_tree(unsigned long ntips=DEFAULT_NUM_TIPS) {
    auto tree = generate_unbalanced_tree(ntips);
    tree.add_coalescent_edge_lengths(0);
    return check_tree(tree, 0, "[UNBALANCED, MEAN-COAL]");
}

int test_converse_coalescent_ages_on_unbalanced_tree(unsigned long ntips=DEFAULT_NUM_TIPS) {
    auto tree = generate_unbalanced_tree(ntips);
    tree.add_coalescent_edge_lengths(1);
    return check_tree(tree, 1, "[UNBALANCED, CONVERSE-COAL]");
}

int test_uniform_coalescent_ages_on_unbalanced_tree(unsigned long ntips=DEFAULT_NUM_TIPS) {
    auto tree = generate_unbalanced_tree(ntips);
    tree.add_coalescent_edge_lengths(2);
    return check_tree(tree, 2, "[UNBALANCED, UNIFORM]");
}

int main(int, const char * argv[]) {
    TEST_DIR = pstrudel::test::get_test_dir(argv[0]);
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "check-coalescent-edge-lengths.py");
    platypus::bind_standard_interface(TREE_WRITER);
    int fails = 0;
    fails += test_mean_coalescent_ages_on_unbalanced_tree();
    fails += test_uniform_coalescent_ages_on_unbalanced_tree();
    fails += test_converse_coalescent_ages_on_unbalanced_tree();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
