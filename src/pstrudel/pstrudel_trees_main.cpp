#include <cassert>
#include <colugo/filesys.hpp>
#include <colugo/cmdopt.hpp>
#include <colugo/logger.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pstrudel.hpp"
#include "distancetree.hpp"

int main(int argc, const char * argv[]) {
    std::string      prog_id                                =     pstrudel::get_program_identification("PSTRUDEL-TREES").c_str();
    std::string      format                                 =     "nexus";
    unsigned         long num_interpolated_points           =     0;
    bool             calculate_baseline_distances           =     false;
    bool             calculate_pairwise_distances           =     false;
    bool             calculate_symmetric_diff               =     false;
    std::string      reference_trees_filepath               =     "";
    std::string      default_output_filename_stem           =     "pstrudel-trees";
    std::string      output_prefix                          =     default_output_filename_stem;
    bool             suppress_copying_of_comparison_trees   =     false;
    bool             suppress_header_row                    =     false;
    bool             replace_existing_output_files          =     false;
    bool             quiet                                  =     false;

    colugo::OptionParser parser = colugo::OptionParser(
            prog_id.c_str(),
            "Calculate and report distances between trees.",
            "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    parser.add_switch(&calculate_baseline_distances,
            "-b",
            "--baseline",
            "Calculate distances from every tree to reference trees"
            " (trees in file given by '--reference-trees' or canonical tree "
            " patterns if '--reference-trees' not specified).");
    parser.add_switch(&calculate_pairwise_distances,
            "-p",
            "--pairwise",
            "Calculate distances from every tree to every other tree.");
    parser.add_switch(&calculate_symmetric_diff,
            "-s",
            "--symmetric-difference",
            "Calculate (unlabeled) Robinson-Foulds symmetric difference distance from"
            " every tree to every other tree.");
    parser.add_option<std::string>(&reference_trees_filepath, "-t", "--reference-trees",
            "Tree(s) to use as benchmark(s) when calculating baseline distances;"
            " if not specified, canonical tree patterns will be used.",
            "REFERENCE-TREE-FILEPATH");
    parser.add_option<std::string>(&format, "-f", "--format",
            "Format for tree sources, one of:'nexus' [default] or 'newick'.",
            "FORMAT");
    parser.add_option<unsigned long>(&num_interpolated_points, "-n", "--profile-size",
            "Number of interpolated points in profile metric; if specified."
            " This must be equal or greater to than the largest input data size."
            " If not specified, will default to the largest input data size.",
            "SIZE");
    parser.add_option<std::string>(&output_prefix, "-o", "--output-prefix",
            "Prefix (directory path and filename stem) for output file(s); if"
            " not specified, will default to '%default'.",
            "PATH/TO/OUTPUT");
    parser.add_switch(&replace_existing_output_files,
            "-r",
            "--replace-existing-output",
            "Replace (overwrite) existing output files.");
    parser.add_switch(&suppress_copying_of_comparison_trees, NULL, "--suppress-comparison-tree-copy",
            "Do not save an aggregated copy of the comparison tree files.");
    parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row",
            "Do not write column/field name row in results.");
    parser.add_switch(&quiet, "-q", "--quiet", "Suppress all informational/progress messages.");
    parser.parse(argc, argv);

    // set up output filepaths
    if (COLUGO_FILESYS_PATH_SEPARATOR[0] == output_prefix[output_prefix.size()-1]) {
        output_prefix = default_output_filename_stem;
    }
    if (output_prefix[output_prefix.size() - 1] != '.') {
        output_prefix = output_prefix + '.';
    }
    std::map<std::string, std::string> output_filepaths;
    output_filepaths["log"] = output_prefix + "log";
    if (!suppress_copying_of_comparison_trees) {
        output_filepaths["comparison-trees"] = output_prefix + "comparison.trees";
    }
    if (calculate_baseline_distances) {
        output_filepaths["reference-trees"] = output_prefix + "reference.trees";
        output_filepaths["baseline-distances"] = output_prefix + "distances.baseline.txt";
        output_filepaths["baseline-distances-stacked"] = output_prefix + "distances.baseline.stacked.txt";
    }
    if (calculate_pairwise_distances) {
        output_filepaths["pairwise-distances"] = output_prefix + "distances.pairwise.txt";
        output_filepaths["pairwise-distances-stacked"] = output_prefix + "distances.pairwise.stacked.txt";
    }

    if (!replace_existing_output_files) {
        // make sure we are not overwriting
        std::set<std::string> existing_output_files;
        for (auto & off : output_filepaths) {
            if (colugo::filesys::exists(off.second)) {
                existing_output_files.insert(off.second);
            }
        }
        if (!existing_output_files.empty()) {
            colugo::console::err_wrapped("The following files will be overwritten by the current command:\n");
            for (auto & ef : existing_output_files) {
                colugo::console::err_line("  - '", ef, "'");
            }
            colugo::console::err_wrapped("\nRe-run the command using the '-r'/'--replace-existing-output' option to overwrite the files "
                    "or specify a different output prefix using the '-o'/'--output-prefix' option. Alternatively, you "
                    "may prefer to rename, move, or simply delete these files before proceeding.");
            exit(EXIT_FAILURE);
        }
    }

    // set up run logger
    colugo::Logger logger("pstrudel-trees");
    if (quiet) {
        logger.add_channel(std::cerr,
                colugo::Logger::LoggingLevel::WARNING,
                true,
                colugo::Logger::LoggingLevel::NOTSET);
    } else {
        logger.add_channel(std::cerr,
                colugo::Logger::LoggingLevel::INFO,
                true,
                colugo::Logger::LoggingLevel::NOTSET);
    }
    std::ofstream logfile(output_filepaths["log"]);
    logger.add_channel(logfile,
            colugo::Logger::LoggingLevel::INFO,
            true,
            colugo::Logger::LoggingLevel::NOTSET);

    // log program identification
    logger.info(prog_id, " is running");

    // get and identify omparison trees
    std::vector<std::string> args = parser.get_args();
    std::vector<pstrudel::DistanceTree> comparison_trees;
    std::vector<std::string> comparison_tree_sources;
    if (args.size() == 0) {
        logger.info("(reading trees for comparison from standard input)");
        pstrudel::treeio::read_from_stream(comparison_trees, std::cin, format);
        comparison_tree_sources.reserve(comparison_trees.size());
        for (unsigned long i = 0; i < comparison_trees.size(); ++i) {
            comparison_tree_sources.push_back("<stdin>");
        }
        logger.info(comparison_trees.size(), " trees read from standard input");
    } else {
        unsigned file_idx = 1;
        unsigned long num_trees_read = 0;
        for (auto & arg : args) {
            auto fullpath = colugo::filesys::absolute_path(arg);
            logger.info("Reading trees for comparison from file ", file_idx, " of ", args.size(), ": '", fullpath, "'");
            num_trees_read = pstrudel::treeio::read_from_filepath(comparison_trees, fullpath, format);
            for (unsigned long i = 0; i < num_trees_read; ++i) {
                comparison_tree_sources.push_back(fullpath);
            }
            logger.info(num_trees_read, " trees read from '", fullpath,
                    "' (total number of trees in comparison set: ", comparison_trees.size(), ")");
            ++file_idx;
        }
        logger.info(comparison_trees.size(), " trees read from ", args.size(), " file(s)");
    }

    // // process source_trees
    // logger.info("Calculating tree metrics ...");
    // for (auto & tree : source_trees) {
    //     tree.calc_profile_metrics();
    // }

    // // get number of interpolated points
    // if (num_interpolated_points == 0) {
    //     logger.info("Calculating number of interpolated points in profiles ...");
    //     std::map<const std::string, unsigned long> global_num_interpolated_points;
    //     for (auto & tree : source_trees) {
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
    //     for (auto & tree : source_trees) {
    //         tree.set_num_interpolated_profile_points(global_num_interpolated_points);
    //     }
    // }

    // // build profiles
    // logger.info("Building tree profiles ...");
    // for (auto & tree : source_trees) {
    //     tree.build_profile();
    // }

    // // calculate pairwise profile distances
    // TreeDistanceMatrixType  unweighted_profile_distances;
    // TreeDistanceMatrixType  weighted_profile_distances;
    // TreeDistanceMatrixType  full_distances;
    // logger.info("Calculating all pairwise profile distances ...");
    // for (unsigned long tree_idx1 = 0; tree_idx1 < source_trees.size() - 1; ++tree_idx1) {
    //     auto & tree1 = source_trees[tree_idx1];
    //     for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < source_trees.size(); ++tree_idx2) {
    //         auto & tree2 = source_trees[tree_idx2];
    //         unweighted_profile_distances[tree_idx1][tree_idx2] = tree1.get_unweighted_subprofile_distance(tree2);
    //         weighted_profile_distances[tree_idx1][tree_idx2] = tree1.get_weighted_subprofile_distance(tree2);
    //         full_distances[tree_idx1][tree_idx2] = tree1.get_distance(tree2);
    //     }
    // }

    // // calculate pairwise symmetric distances
    // TreeDistanceMatrixType  unlabeled_symmetric_differences;
    // if (calculate_other_distance_metrics) {
    //     logger.info("Calculating other metrics ...");
    //     std::vector<pstrudel::SymmetricDifferenceTree> symmetric_difference_trees;
    //     symmetric_difference_trees.reserve(source_trees.size());
    //     for (auto & tree1 : source_trees) {
    //         symmetric_difference_trees.emplace_back(tree1);
    //         auto & sd_tree = symmetric_difference_trees.back();
    //         sd_tree.calc_subtree_sizes();
    //     }
    //     COLUGO_ASSERT(symmetric_difference_trees.size() == source_trees.size());
    //     for (unsigned long tree_idx1 = 0; tree_idx1 < symmetric_difference_trees.size() - 1; ++tree_idx1) {
    //         auto & tree1 = symmetric_difference_trees[tree_idx1];
    //         for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < symmetric_difference_trees.size(); ++tree_idx2) {
    //             auto & tree2 = symmetric_difference_trees[tree_idx2];
    //             unlabeled_symmetric_differences[tree_idx1][tree_idx2] = tree1.calc_leaf_set_sizes_unlabeled_symmetric_difference(tree2);
    //         }
    //     }
    // }

    // std::ostream& out = std::cout;
    // if (!suppress_header_row) {
    //     out << "Tree_i";
    //     out << "\t" << "Tree_j";
    //     out << "\t" << "uw.profile.dist";
    //     out << "\t" << "w.profile.dist";
    //     out << "\t" << "a.profile.dist";
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << "sym.diff.dist";
    //     }
    //     out << std::endl;
    // }

    // unsigned long num_trees = source_trees.size();
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


