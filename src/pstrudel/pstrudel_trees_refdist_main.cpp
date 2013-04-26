#include <cassert>
#include <colugo/cmdopt.hpp>
#include <colugo/utility.hpp>
#include <colugo/logger.hpp>
#include <colugo/textutil.hpp>
#include <platypus/model/archetypaltree.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pairwise_distance_tree.hpp"
#include "symmetric_difference_tree.hpp"
#include "pstrudel.hpp"

typedef std::map<unsigned long, double> TreeReferenceDistanceType;

int main(int argc, const char * argv[]) {
    std::string prog_id = pstrudel::get_program_identification("PSTRUDEL-TREES-REFDIST").c_str();
    std::string format = "nexus";
    // bool force_rooted = false;
    // bool force_unrooted = false;
    unsigned long num_interpolated_points = 0;
    bool calculate_other_distance_metrics = false;
    std::string output_prefix = "";
    bool suppress_header_row = false;
    bool quiet = false;

    colugo::OptionParser parser = colugo::OptionParser(
            prog_id.c_str(),
            "Calculate and report distances between each tree and reference trees (maximally-balanced and maximally-unbalanced trees).",
            "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    parser.add_option<std::string>(&format, "-f", "--format",
            "format for input source ('nexus', 'newick'; default = 'nexus')");
    parser.add_option<unsigned long>(&num_interpolated_points, "-n", "--profile-size",
            "number of interpolated points in profile; if specified, this must be equal or greater to than the largest input data size; if not specified, will default to the largest input data size");
    parser.add_switch(&calculate_other_distance_metrics,
            "-a",
            "--calculate-other-distance-metrics",
            "calculate other distance metrics (i.e., unlabeled symmetric difference, etc.)");
    parser.add_option<std::string>(&output_prefix, "-o", "--output-prefix",
            "prefix (directory path and filename stem) for output file(s); if not specified, will write to standard output");
    parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row", "do not write column/field name row in reuslts");
    // parser.add_switch(&rooted, "--rooted", "force all pairwise_distance_trees to be rooted");
    // parser.add_switch(&unrooted, "--unrooted", "force all pairwise_distance_trees to be unrooted");
    parser.add_switch(&quiet, "-q", "--quiet", "suppress all informational/progress messages");
    parser.parse(argc, argv);

    // set up logger
    colugo::Logger logger("pstrudel-trees-refdist");
    if (quiet) {
        logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::WARNING);
    } else {
        logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::INFO);
    }

    // get source pairwise_distance_trees
    std::vector<std::string> args = parser.get_args();
    std::vector<pstrudel::PairwiseDistanceTree> pairwise_distance_trees;
    if (args.size() == 0) {
        logger.info("(reading trees from standard input)");
        pstrudel::treeio::read_from_stream(pairwise_distance_trees, std::cin, format);
    } else {
        unsigned file_idx = 1;
        for (auto & arg : args) {
            logger.info("Reading tree source file ", file_idx, " of ", args.size(), ": ", arg);
            pstrudel::treeio::read_from_filepath(pairwise_distance_trees, arg, format);
            ++file_idx;
        }
        logger.info(pairwise_distance_trees.size(), " trees read.");
    }

    // process pairwise_distance_trees
    logger.info("Calculating tree metrics ...");
    for (auto & tree : pairwise_distance_trees) {
        tree.calc_profile_metrics();
    }

    // get number of interpolated points
    std::map<const std::string, unsigned long> global_num_interpolated_points;
    if (num_interpolated_points == 0) {
        logger.info("Calculating number of interpolated points in profiles ...");
        for (auto & tree : pairwise_distance_trees) {
            tree.poll_max_num_interpolated_profile_points(global_num_interpolated_points);
        }
        unsigned long min_points = 0;
        unsigned long max_points = 0;
        for (auto npi : global_num_interpolated_points) {
            if (min_points < npi.second) {
                min_points = npi.second;
            }
            if (max_points < npi.second) {
                max_points = npi.second;
            }
        }
        logger.info("Minumum number of interpolated profile points: ", min_points);
        logger.info("Maximum number of interpolated profile points: ", max_points);
        logger.info("Configuring profile build regime ...");
        for (auto & tree : pairwise_distance_trees) {
            tree.set_num_interpolated_profile_points(global_num_interpolated_points);
        }
    }

    // build profiles
    logger.info("Building tree profiles ...");
    for (auto & tree : pairwise_distance_trees) {
        tree.build_profile();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // TODO: handle this more smartly
    // for now, assume all trees have the same number of leaves
    // generation of reference trees
    unsigned long num_leaves = pairwise_distance_trees[0].get_num_leaves();
    std::vector<pstrudel::PairwiseDistanceNodeValue> leaves;

    // balanced tree
    logger.info("Generating balanced tree ..");
    leaves.reserve(num_leaves);
    for (unsigned long i = 0; i < num_leaves; ++i) {
        leaves.emplace_back();
    }
    pstrudel::PairwiseDistanceTree balanced_tree;
    build_maximally_balanced_tree(balanced_tree, leaves.begin(), leaves.end());
    // std::ostringstream balanced_tree_newick;
    // pstrudel::treeio::write_nexus(balanced_tree, balanced_tree_newick, false);
    // logger.info("Balanced tree: ", balanced_tree_newick.str());
    balanced_tree.calc_profile_metrics();
    if (num_interpolated_points == 0) {
        balanced_tree.set_num_interpolated_profile_points(global_num_interpolated_points);
    }

    leaves.clear();
    logger.info("Generating unbalanced tree ..");
    for (unsigned long i = 0; i < num_leaves; ++i) {
        leaves.emplace_back();
    }
    pstrudel::PairwiseDistanceTree unbalanced_tree;
    build_maximally_unbalanced_tree(unbalanced_tree, leaves.begin(), leaves.end());
    // std::ostringstream unbalanced_tree_newick;
    // pstrudel::treeio::write_nexus(unbalanced_tree, unbalanced_tree_newick, false);
    // logger.info("Balanced tree: ", unbalanced_tree_newick.str());
    unbalanced_tree.calc_profile_metrics();
    if (num_interpolated_points == 0) {
        unbalanced_tree.set_num_interpolated_profile_points(global_num_interpolated_points);
    }

    // calculate pairwise profile distances
    TreeReferenceDistanceType  balanced_tree_profile_distances;
    TreeReferenceDistanceType  unbalanced_tree_profile_distances;
    logger.info("Calculating all pairwise profile distances ...");
    for (unsigned long tree_idx1 = 0; tree_idx1 < pairwise_distance_trees.size(); ++tree_idx1) {
        auto & tree1 = pairwise_distance_trees[tree_idx1];
        balanced_tree_profile_distances[tree_idx1] = tree1.get_unweighted_subprofile_distance(balanced_tree);
        unbalanced_tree_profile_distances[tree_idx1] = tree1.get_unweighted_subprofile_distance(unbalanced_tree);
        // for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < pairwise_distance_trees.size(); ++tree_idx2) {
        //     auto & tree2 = pairwise_distance_trees[tree_idx2];
        //     unweighted_profile_distances[tree_idx1][tree_idx2] = tree1.get_unweighted_subprofile_distance(tree2);
        //     weighted_profile_distances[tree_idx1][tree_idx2] = tree1.get_weighted_subprofile_distance(tree2);
        //     full_distances[tree_idx1][tree_idx2] = tree1.get_distance(tree2);
        // }
    }

    // calculate pairwise symmetric distances
    // TreeReferenceDistanceType  balanced_tree_symdiff_distances;
    // TreeReferenceDistanceType  unbalanced_tree_symdiff_distances;
    // if (calculate_other_distance_metrics) {
    //     logger.info("Calculating other metrics ...");
    //     std::vector<pstrudel::SymmetricDifferenceTree> symmetric_difference_trees;
    //     symmetric_difference_trees.reserve(pairwise_distance_trees.size());
    //     for (auto & tree1 : pairwise_distance_trees) {
    //         symmetric_difference_trees.emplace_back(tree1);
    //         auto & sd_tree = symmetric_difference_trees.back();
    //         sd_tree.calc_subtree_sizes();
    //     }
    //     COLUGO_ASSERT(symmetric_difference_trees.size() == pairwise_distance_trees.size());
    //     for (unsigned long tree_idx1 = 0; tree_idx1 < symmetric_difference_trees.size(); ++tree_idx1) {
    //         auto & tree1 = symmetric_difference_trees[tree_idx1];
    //         balanced_tree_symdiff_distances[tree_idx1] = calc_leaf_set_sizes_unlabeled_symmetric_difference(balanced_tree2);
    //         unbalanced_tree_symdiff_distances[tree_idx1] = calc_leaf_set_sizes_unlabeled_symmetric_difference(unbalanced_tree2);
    //         // for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < symmetric_difference_trees.size(); ++tree_idx2) {
    //         //     auto & tree2 = symmetric_difference_trees[tree_idx2];
    //         //     unlabeled_symmetric_differences[tree_idx1][tree_idx2] = tree1.calc_leaf_set_sizes_unlabeled_symmetric_difference(tree2);
    //         // }
    //     }
    // }

    // std::ostream& out = std::cout;
    // if (!suppress_header_row) {
    //     out << "Tree_i";
    //     out << "\t" << "Tree_j";
    //     out << "\t" << "Unweighted.Profile.Distance";
    //     out << "\t" << "Edge.Weighted.Profile.Distance";
    //     out << "\t" << "Aggregate.Profile.Distance";
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << "Unlabeled.Symmetric.Difference.Sets";
    //     }
    //     out << std::endl;
    // }

    // unsigned long num_trees = pairwise_distance_trees.size();
    // for (unsigned long tree_idx1 = 0; tree_idx1 < num_trees - 1; ++tree_idx1) {
    //     for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < num_trees; ++tree_idx2) {
    //         out << tree_idx1 + 1;
    //         out << "\t" << tree_idx2 + 1;
    //         out << "\t" << std::setprecision(20) << unweighted_profile_distances[tree_idx1][tree_idx2];
    //         out << "\t" << std::setprecision(20) << weighted_profile_distances[tree_idx1][tree_idx2];
    //         out << "\t" << std::setprecision(20) << full_distances[tree_idx1][tree_idx2];
    //         if (calculate_other_distance_metrics) {
    //             out << "\t" << unlabeled_symmetric_differences[tree_idx1][tree_idx2];
    //         }
    //         out << std::endl;
    //     }
    // }

}


