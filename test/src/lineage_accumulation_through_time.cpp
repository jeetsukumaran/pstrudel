#include <vector>
#include <iostream>
#include <sstream>
#include <colugo/subprocess.hpp>
#include <platypus/serialize/newick.hpp>
#include <pstrudel/dataio.hpp>
#include <pstrudel/distancetree.hpp>
#include "pstrudel_testing.hpp"

std::string TEST_DIR;
std::string CHECK_SCRIPT;
static const unsigned long DEFAULT_NUM_TIPS = 50;
platypus::NewickWriter<pstrudel::DistanceTree> TREE_WRITER;

int test_file(const std::string & test_data_filepath,
        const std::string & label_prefix,
        bool verify=true) {
    int fails = 0;
    std::vector<pstrudel::DistanceTree>  trees;
    pstrudel::treeio::read_from_filepath(trees, test_data_filepath, "nexus");
    unsigned long tree_idx = 0;
    for (auto & tree : trees) {
        std::string label = label_prefix + ":" + std::to_string(tree_idx + 1);
        pstrudel::LineageThroughTimeProfileCalculator ltt_calc(tree);
        auto num_transects = ltt_calc.get_default_num_transects();
        auto transect_offsets = ltt_calc.build_transect_offsets(num_transects);
        assert(transect_offsets.size() > 0);
        auto ltt_profile = ltt_calc.build_lineage_accumulation_through_time_profile(transect_offsets);
        if (verify) {
            std::ostringstream o;
            TREE_WRITER.write(o, tree);
            o << "\n";
            for (unsigned long idx = 0; idx < transect_offsets.size(); ++idx) {
                o << idx << "\t";
                o << std::fixed << std::setprecision(22) << transect_offsets[idx] << "\t";
                o << std::fixed << std::setprecision(22) << ltt_profile.raw_data(idx);
                o << "\n";
            }
            // std::cerr << label << " >>>" << "\n";
            // std::cerr << o.str() << std::endl;
            colugo::Subprocess ps({"python", CHECK_SCRIPT, "-f", "newick", "-c", std::to_string(num_transects), "-l", label});
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
                // std::string out;
                // std::string err;
                // ps.read_pipes(out, err);
                // std::cerr << out;
                // std::cerr << err;
                // std::cerr << "---" << std::endl;
                // std::cout << o.str() << std::endl;
                // std::cerr << "---" << std::endl;
                // std::cerr << "-f newick -c " << std::to_string(num_transects) << std::endl;
                // std::cerr << "---" << std::endl;
                exit(1);
            }
        } else {
            // don't verify: called when we are testing tree w/o branch lengths
            // just to make sure that the calculations proceed OK
        }
        ++tree_idx;
    }
    return fails;
}

int test_ltt1() {
    std::string file_basename =  "pythonidae.reference-trees.nexus";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix, true);
}

int test_ltt2() {
    std::string file_basename =  "apternodus.tre";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix, false);
}

int test_ltt3() {
    std::string file_basename =  "pythonidae.mb.run1.t";
    std::string test_data_filepath = pstrudel::test::join_path(TEST_DIR, "data", "trees", "general", file_basename);
    std::string label_prefix = file_basename;
    return test_file(test_data_filepath, label_prefix, false);
}

int main(int, const char * argv[]) {
    TEST_DIR = pstrudel::test::get_test_dir(argv[0]);
    CHECK_SCRIPT = pstrudel::test::join_path(TEST_DIR, "check-lineage-accumulation-through-time.py");
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
