#include <vector>
#include <string>
#include <iostream>
#include <colugo/subprocess.hpp>
#include <platypus/serialize/newick.hpp>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;
platypus::NewickWriter<pstrudel::DistanceTree> TREE_WRITER;

int test_file(const std::string & test_data_filename,
        const std::string & format,
        const std::string label_prefix) {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", test_data_filename);
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, format);
    unsigned long tree_idx1 = 0;
    unsigned long fails = 0;
    for (auto & tree : trees) {
        std::string label = label_prefix + ":" + std::to_string(tree_idx1 + 1);
        std::ostringstream o;
        TREE_WRITER.write(o, tree);
        o << "\n";

        // check node ages
        auto node_ages = tree.calc_node_ages(false);
        for (auto & na : node_ages) {
            o << na << ",";
        }
        o << "\n";

        o << std::fixed << std::setprecision(22) << tree.get_pybus_harvey_gamma();
        o << "\n";
        colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick", "-l", label});
        auto result = ps.communicate(o.str());
        if (ps.returncode() != 0) {
            std::cerr << "\n-- test '" << label << "' returned error: " << ps.returncode() << "\n";
            std::cerr << result.first;
            std::cerr << result.second;
            fails += 1;
        }
        ++tree_idx1;
    }
    return fails;
}

int test_dist1() {
    int fails = 0;
    std::string src_filepath = "n10-unbalanced.nexus.trees";
    return test_file(src_filepath, "nexus", src_filepath);
}

int test_dist2() {
    int fails = 0;
    std::string src_filepath = "pythonidae.reference-trees.nexus";
    return test_file(src_filepath, "nexus", src_filepath);
}

int main(int, const char * argv[]) {
    TEST_DIR = pstrudel::test::get_test_dir(argv[0]);
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "calc-tree-unary-stats.py");
    platypus::bind_standard_interface(TREE_WRITER);
    TREE_WRITER.set_edge_length_precision(22);
    int fails = 0;
    fails += test_dist1();
    // fails += test_dist2();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
