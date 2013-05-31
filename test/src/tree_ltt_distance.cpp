#include <vector>
#include <string>
#include <iostream>
#include <colugo/subprocess.hpp>
#include <platypus/serialize/newick.hpp>
#include <pstrudel/treeshape.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;
platypus::NewickWriter<pstrudel::TreeShape> TREE_WRITER;

int test_file(const std::string & test_data_filename,
        const std::string & format,
        const std::string label_prefix) {
    typedef pstrudel::TreeShape  TreeType;
    std::vector<TreeType>  trees;
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", test_data_filename);
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, format);
    unsigned long tree_idx1 = 0;
    unsigned long tree_idx2 = 0;
    unsigned long fails = 0;
    for (auto & tree1 : trees) {
        for (auto & tree2 : trees) {
            std::string label = label_prefix + ":" + std::to_string(tree_idx1 + 1) + ":" + std::to_string(tree_idx2 + 1);
            std::ostringstream o;
            TREE_WRITER.write(o, tree1);
            o << "\n";
            TREE_WRITER.write(o, tree2);
            o << "\n";
            o << std::fixed << std::setprecision(22) << tree1.get_lineage_accumulation_profile_distance(tree2);
            o << "\n";
            o << std::fixed << std::setprecision(22) << tree1.get_lineage_splitting_time_profile_distance(tree2);
            o << "\n";
            o << std::fixed << std::setprecision(22) << tree1.get_scaled_lineage_splitting_time_profile_distance(tree2);
            o << "\n";
            std::vector<pstrudel::Profile *> profiles;
            auto & calc1 = tree1.get_lineage_through_time_calculator();
            auto & calc2 = tree2.get_lineage_through_time_calculator();
            for (auto & calc : {&calc1, &calc2}) {
                profiles.push_back(&(calc->get_lineage_accumulation_through_time_profile()));
                profiles.push_back(&(calc->get_lineage_splitting_time_profile()));
                profiles.push_back(&(calc->get_scaled_lineage_splitting_time_profile()));
            }
            for (auto & profile : profiles) {
                for (auto & val : profile->get_profile()) {
                    o << val << ",";
                }
                o << "\n";
            }
            colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick", "-l", label});
            auto result = ps.communicate(o.str());
            if (ps.returncode() != 0) {
                std::cerr << "\n-- test '" << label << "' returned error: " << ps.returncode() << "\n";
                std::cerr << result.first;
                std::cerr << result.second;
                fails += 1;
            }
            ++tree_idx2;
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
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "calc-tree-ltt-distance.py");
    platypus::bind_standard_interface(TREE_WRITER);
    TREE_WRITER.set_edge_length_precision(22);
    int fails = 0;
    fails += test_dist1();
    fails += test_dist2();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
