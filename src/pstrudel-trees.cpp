#include "colugo-utilities/src/cmdopt.hpp"
#include "colugo-utilities/src/utility.hpp"
#include "colugo-utilities/src/logger.hpp"
#include "colugo-utilities/src/textutil.hpp"
#include "dataio.hpp"
#include "split.hpp"
#include "pairwise_distance_tree.hpp"
#include "pstrudel.hpp"

int main(int argc, const char * argv[]) {
    // std::string prog_id = pstrudel::get_program_identification("PSTRUDEL-TREES").c_str();
    // std::string format = "nexus";
    // // bool force_rooted = false;
    // // bool force_unrooted = false;
    // unsigned long num_interpolated_points = 0;
    // bool calculate_other_distance_metrics = false;
    // // bool calculate_tree_metrics = false;
    // std::string output_prefix = "";
    // bool suppress_header_row = false;
    // bool quiet = false;

    // colugo::OptionParser parser = colugo::OptionParser(
    //         prog_id.c_str(),
    //         "Calculate and report distances between every distinct pair of trees.",
    //         "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    // parser.add_option<std::string>(&format, "-f", "--format",
    //         "format for input source ('nexus', 'newick'; default = 'nexus')");
    // parser.add_option<unsigned long>(&num_interpolated_points, "-n", "--profile-size",
    //         "number of interpolated points in profile; if specified, this must be equal or greater to than the largest input data size; if not specified, will default to the largest input data size");
    // // parser.add_switch(&calculate_other_distance_metrics,
    // //         "-s",
    // //         "--report-tree-metrics",
    // //         "report tree metrics (e.g., number of tips, sum of tree edge lengths)");
    // parser.add_switch(&calculate_other_distance_metrics,
    //         "-a",
    //         "--calculate-other-distance-metrics",
    //         "calculate other distance metrics (e.g., symmetric difference, weighted Robinson-Foulds distance, edge length Euclidean distances)");
    // parser.add_option<std::string>(&output_prefix, "-o", "--output-prefix",
    //         "prefix (directory path and filename stem) for output file(s); if not specified, will write to standard output");
    // parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row", "do not write column/field name row in reuslts");
    // // parser.add_switch(&rooted, "--rooted", "force all trees to be rooted");
    // // parser.add_switch(&unrooted, "--unrooted", "force all trees to be unrooted");
    // parser.add_switch(&quiet, "-q", "--quiet", "suppress all informational/progress messages");
    // parser.parse(argc, argv);

    // // set up logger
    // colugo::Logger logger("pstrudel-trees");
    // if (quiet) {
    //     logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::WARNING);
    // } else {
    //     logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::INFO);
    // }

    // // get source trees
    // std::vector<std::string> args = parser.get_args();
    // std::vector<pstrudel::PairwiseDistanceTree> trees;
    // if (args.size() == 0) {
    //     logger.info("(reading trees from standard input)");
    //     pstrudel::treeio::read_from_stream(trees, std::cin, format);
    // } else {
    //     unsigned file_idx = 1;
    //     for (auto & arg : args) {
    //         logger.info("Reading tree source file ", file_idx, " of ", args.size(), ": ", arg);
    //         pstrudel::treeio::read_from_filepath(trees, arg, format);
    //         ++file_idx;
    //     }
    //     logger.info(trees.size(), " trees read.");
    // }

    // // process trees
    // logger.info("Calculating tree metrics and building tree profiles ...");
    // pstrudel::TreeComparisonCalculator tree_comparison_calculator;
    // for (auto & tree : trees) {
    //     tree.calc_profile_metrics();
    //     if (calculate_other_distance_metrics) {
    //         tree_comparison_calculator.calc_tree_splits(tree);
    //     }
    // }

    // // get number of interpolated points
    // if (num_interpolated_points == 0) {
    //     logger.info("Calculating number of interpolated points in profiles ...");
    //     std::map<const std::string, unsigned long> global_num_interpolated_points;
    //     for (auto & tree : trees) {
    //         tree.poll_max_num_interpolated_profile_points(global_num_interpolated_points);
    //     }
    //     unsigned long min_points = 0;
    //     unsigned long max_points = 0;
    //     for (auto npi : global_num_interpolated_points) {
    //         if (min_points < npi.second) {
    //             min_points = npi.second;
    //         }
    //         if (max_points < npi.second) {
    //             max_points = npi.second;
    //         }
    //     }
    //     logger.info("Minumum number of interpolated profile points: ", min_points);
    //     logger.info("Maximum number of interpolated profile points: ", max_points);
    //     logger.info("Configuring profile build regime ...");
    //     for (auto & tree : trees) {
    //         tree.set_num_interpolated_profile_points(global_num_interpolated_points);
    //     }
    // }

    // // build profiles
    // logger.info("Building tree profiles ...");
    // for (auto & tree : trees) {
    //     tree.build_profile();
    // }

    // std::ostream& out = std::cout;
    // if (!suppress_header_row) {
    //     out << "Tree_i";
    //     out << "\t" << "Tree_j";
    //     out << "\t" << "Unweighted.Profile.Distance";
    //     out << "\t" << "Edge.Weighted.Profile.Distance";
    //     out << "\t" << "Aggregate.Profile.Distance";
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << "Symmetric.Difference";
    //         out << "\t" << "Edge.Absolute.Distance";
    //         out << "\t" << "Edge.Euclidean.Distance";
    //     }
    //     out << std::endl;
    // }

    // for (unsigned long tree_idx1 = 0; tree_idx1 < trees.size() - 1; ++tree_idx1) {
    //     auto & tree1 = trees[tree_idx1];
    //     for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < trees.size(); ++tree_idx2) {
    //         auto & tree2 = trees[tree_idx2];
    //         out << tree_idx1 + 1;
    //         out << "\t" << tree_idx2 + 1;
    //         out << "\t" << std::setprecision(20) << tree1.get_unweighted_subprofile_distance(tree2);
    //         out << "\t" << std::setprecision(20) << tree1.get_weighted_subprofile_distance(tree2);
    //         out << "\t" << std::setprecision(20) << tree1.get_distance(tree2);
    //         if (calculate_other_distance_metrics) {
    //             auto result = tree_comparison_calculator.compare_trees(tree1, tree2);
    //             std::cout << "\t" << result.symmetric_difference;
    //             std::cout << "\t" << std::setprecision(20) << result.edge_len_euclidean_distance;
    //             std::cout << "\t" << std::setprecision(20) << result.weighted_rf_distance;
    //         }
    //         out << std::endl;
    //     }
    // }

}


