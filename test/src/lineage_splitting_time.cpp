#include <vector>
#include <iostream>
#include <sstream>
#include <colugo/subprocess.hpp>
#include <platypus/serialize/newick.hpp>
#include <pstrudel/dataio.hpp>
#include <pstrudel/treeshape.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;
static const unsigned long DEFAULT_NUM_TIPS = 50;
platypus::NewickWriter<pstrudel::TreeShape> TREE_WRITER;

int test_file(const std::string & test_data_filepath,
        const std::string & label_prefix) {
    int fails = 0;
    std::vector<pstrudel::TreeShape>  trees;
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");
    unsigned long tree_idx = 0;
    for (auto & tree : trees) {
        std::string label = label_prefix + ":" + std::to_string(tree_idx + 1);
        pstrudel::LineageThroughTimeProfileCalculator ltt_calc(tree);
        auto profiles = ltt_calc.build_lineage_splitting_time_profile();
        auto lst_profile = profiles.first;
        auto scaled_lst_profile = profiles.second;
        std::ostringstream o;
        TREE_WRITER.write(o, tree);
        o << "\n";
        assert(lst_profile.data_size() == scaled_lst_profile.data_size());
        auto rds = lst_profile.data_size();
        for (unsigned long idx = 0; idx < rds; ++idx) {
            o << std::fixed << std::setprecision(22) << lst_profile.raw_data(idx);
            o << "\t" << std::fixed << std::setprecision(22) << scaled_lst_profile.raw_data(idx);
            o << "\n";
        }
        colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick", "-l", label});
        try {
            auto result = ps.communicate(o.str());
            if (ps.returncode() != 0) {
                std::cerr << "(test '" << label << "' returned error: " << ps.returncode() << ")\n";
                // TREE_WRITER.write(std::cerr, tree);
                // std::cerr << std::endl;
                std::cerr << result.first;
                std::cerr << result.second;
                fails += 1;
            // } else {
            //     std::cerr << result.second;
            }
        } catch (const colugo::SubprocessTimeOutException & e) {
            std::cerr << "(test '" << label << "' timed out)\n";
            exit(1);
        }
        ++tree_idx;
    }
    return fails;
}

int test_ltt1() {
    std::string file_basename =  "pythonidae.reference-trees.nexus";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix);
}

int test_ltt2() {
    std::string file_basename =  "apternodus.tre";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix);
}

int test_ltt3() {
    std::string file_basename =  "pythonidae.mb.run1.t";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix);
}

int main(int, const char * argv[]) {
    TEST_DIR = pstrudel::test::get_test_dir(argv[0]);
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "check-lineage-spltting-times.py");
    platypus::bind_standard_interface(TREE_WRITER);
    TREE_WRITER.set_edge_length_precision(22);
    int fails = 0;
    fails += test_ltt1();
    fails += test_ltt2();
    fails += test_ltt3();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
