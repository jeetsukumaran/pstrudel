#include <vector>
#include <string>
#include <iostream>
#include <colugo/subprocess.hpp>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;

int test_file(const std::string & test_data_filename, const std::string & format, const std::string label) {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", test_data_filename);
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, format);

    std::map<unsigned long, std::map<unsigned long, double>> lc_dists;
    std::map<unsigned long, std::map<unsigned long, double>> lst_dists;
    std::map<unsigned long, std::map<unsigned long, double>> slst_dists;
    std::set<std::pair<unsigned long, unsigned long>> comparisons;
    for (unsigned tidx1 = 0; tidx1 < trees.size(); ++tidx1) {
        for (unsigned tidx2 = 0; tidx2 < trees.size(); ++tidx2) {
            auto & tree1 = trees[tidx1];
            auto & tree2 = trees[tidx2];
            lc_dists[tidx1][tidx2] = tree1.get_lineage_accumulation_profile_distance(tree2);
            lst_dists[tidx1][tidx2] = tree1.get_lineage_splitting_time_profile_distance(tree2);
            slst_dists[tidx1][tidx2] = tree1.get_scaled_lineage_splitting_time_profile_distance(tree2);
            comparisons.insert(std::make_pair(tidx1, tidx2));
        }
    }

    colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", format, "-l", label, test_data_filepath, "--verbosity", "10"});
    int retcode = ps.wait();
    auto ps_stdout = ps.get_stdout();
    auto ps_stderr = ps.get_stderr();
    if (retcode != 0) {
        std::cerr << "(test '" << label << "' returned error: " << retcode << ")\n";
        std::cerr << ps_stdout;
        std::cerr << ps_stderr;
        return 1;
    }
    int fails = 0;
    unsigned long tidx1   = 0;
    unsigned long tidx2   = 0;
    double exp_lc_d   = 0.0;
    double exp_lst_d  = 0.0;
    double exp_slst_d = 0.0;
    std::istringstream src(ps_stdout);
    while (src.good()) {
        src >> tidx1;
        src >> tidx2;
        src >> exp_lc_d;
        src >> exp_lst_d;
        src >> exp_slst_d;
        fails += platypus::testing::compare_almost_equal(
                exp_lc_d,
                lc_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect lineage accumulation through time profile distance for tree ", tidx1, " and tree ", tidx2);
        fails += platypus::testing::compare_almost_equal(
                exp_lst_d,
                lst_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect lineage splitting times profile distance for tree ", tidx1, " and tree ", tidx2);
        fails += platypus::testing::compare_almost_equal(
                exp_slst_d,
                slst_dists[tidx1][tidx2],
                __FILE__,
                __LINE__,
                "Incorrect scaled lineage splitting times profile distance for tree ", tidx1, " and tree ", tidx2);
        comparisons.erase(std::make_pair(tidx1, tidx2));
    }
    fails += platypus::testing::compare_equal(
            0UL,
            comparisons.size(),
            __FILE__,
            __LINE__,
            "Not all trees checked");
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
    int fails = 0;
    fails += test_dist1();
    fails += test_dist2();
    if (fails != 0) {
        return EXIT_FAILURE;
    } else {
        return EXIT_SUCCESS;
    }
}
