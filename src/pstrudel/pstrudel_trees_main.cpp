#include <cassert>
#include <colugo/filesys.hpp>
#include <colugo/cmdopt.hpp>
#include <colugo/logger.hpp>
#include <platypus/model/datatable.hpp>
#include <platypus/serialize/newick.hpp>
#include <platypus/model/treepattern.hpp>
#include <platypus/model/standardinterface.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pstrudel.hpp"
#include "distancetree.hpp"

typedef std::map<unsigned long, std::map<std::string, pstrudel::DistanceTree>> TreePatternCollectionType;
typedef std::map<unsigned long, std::map<std::string, double>> TreePatternMaxDistances;

template <class TreeT>
std::vector<typename TreeT::value_type> generate_leaves(unsigned long num_tips,
        bool labeled) {
    std::vector<typename TreeT::value_type> leaves;
    leaves.reserve(num_tips);
    for (unsigned long i = 0; i < num_tips; ++i) {
        leaves.emplace_back();
        if (labeled) {
            auto & leaf = leaves.back();
            leaf.set_label("t" + std::to_string(i+1));
        }
    }
    return leaves;
}

template <typename BuildFnT>
void build_tree_pattern(pstrudel::DistanceTree & tree, BuildFnT tree_building_fn, unsigned long num_tips, bool labeled) {
    auto leaves = generate_leaves<pstrudel::DistanceTree>(num_tips, labeled);
    tree_building_fn(tree, leaves.begin(), leaves.end());
}

void build_canonical_tree_patterns(
        unsigned long num_tips,
        TreePatternCollectionType & tree_patterns,
        TreePatternMaxDistances & tree_pattern_max_unweighted_pairwise_tip_profile_distances,
        TreePatternMaxDistances & tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances,
        bool labeled) {
    build_tree_pattern(
            tree_patterns[num_tips]["max.unbalanced"],
            [](pstrudel::DistanceTree & tree,
                std::vector<pstrudel::DistanceNodeValue>::iterator i1,
                std::vector<pstrudel::DistanceNodeValue>::iterator i2){ platypus::build_maximally_unbalanced_tree(tree, i1, i2);},
            num_tips,
            labeled);
    build_tree_pattern(
            tree_patterns[num_tips]["max.balanced"],
            [](pstrudel::DistanceTree & tree,
                std::vector<pstrudel::DistanceNodeValue>::iterator i1,
                std::vector<pstrudel::DistanceNodeValue>::iterator i2){ platypus::build_maximally_balanced_tree(tree, i1, i2);},
            num_tips,
            labeled);
    double max_dist1 = tree_patterns[num_tips]["max.unbalanced"].get_unweighted_pairwise_tip_profile_distance(tree_patterns[num_tips]["max.balanced"]);
    tree_pattern_max_unweighted_pairwise_tip_profile_distances[num_tips]["max.unbalanced"] = max_dist1;
    tree_pattern_max_unweighted_pairwise_tip_profile_distances[num_tips]["max.balanced"] = max_dist1;
    double max_dist2 = tree_patterns[num_tips]["max.unbalanced"].get_unlabeled_symmetric_difference(tree_patterns[num_tips]["max.balanced"]);
    tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances[num_tips]["max.unbalanced"] = max_dist2;
    tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances[num_tips]["max.balanced"] = max_dist2;
    // auto & ub_tree = tree_patterns[num_tips]["max.unbalanced"];
    // auto ub_leaves = generate_leaves<pstrudel::DistanceTree>(num_tips, labeled);
    // platypus::build_maximally_unbalanced_tree(ub_tree, ub_leaves.begin(), ub_leaves.end());
    // ub_tree.build_pairwise_tip_distance_profiles();
    // ub_tree.calc_subtree_sizes();
    // auto & bal_tree = tree_patterns[num_tips]["max.balanced"];
    // auto bal_leaves = generate_leaves<pstrudel::DistanceTree>(num_tips, labeled);
    // platypus::build_maximally_balanced_tree(bal_tree, bal_leaves.begin(), bal_leaves.end());
    // bal_tree.build_pairwise_tip_distance_profiles();
    // bal_tree.calc_subtree_sizes();
}

int main(int argc, const char * argv[]) {
    std::string      prog_name                               = "pstrudel-trees";
    std::string      prog_id                                 = pstrudel::get_program_identification(prog_name).c_str();
    std::string      format                                  = "nexus";
    unsigned         long num_interpolated_points            = 0;
    bool             calculate_reference_distances           = false;
    bool             calculate_pairwise_distances            = false;
    bool             calculate_symmetric_diff                = false;
    std::string      reference_trees_filepath                = "";
    std::string      default_output_filename_stem            = "pstrudel-trees-results";
    std::string      output_prefix                           = default_output_filename_stem;
    bool             create_aggregated_comparison_trees_copy = false;
    bool             suppress_header_row                     = false;
    bool             add_tree_source_key                     = false;
    bool             replace_existing_output_files           = false;
    unsigned long    log_frequency                           = 0;
    bool             quiet                                   = false;

    colugo::OptionParser parser = colugo::OptionParser(
            prog_id.c_str(),
            "Calculate structural distances between unlabeled phylogenetic trees of various sizes.",
            "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    parser.add_switch(&calculate_reference_distances,
            "-r",
            "--reference",
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
            " every tree to every other tree (within the comparison set, to the"
            " reference trees, or both, depending on whether '-r'/'--reference', "
            " '-p'/'--pairwise', or both options are specified).");
    parser.add_option<std::string>(&reference_trees_filepath, "-t", "--reference-trees",
            "Tree(s) to use as benchmark(s) when calculating reference distances;"
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
            "-x",
            "--replace-existing-output",
            "Replace (overwrite) existing output files.");
    parser.add_switch(&create_aggregated_comparison_trees_copy, NULL, "--save-comparison-trees",
            "Save a copy of the comparison trees in a file (aggregating them from across multiple files if multiple source files specified).");
    parser.add_switch(&add_tree_source_key, NULL, "--add-tree-source-key",
            "Add a column in the results identifying the filename of the source of the tree(s) being compared.");
    parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row",
            "Do not write column/field name row in results.");
    parser.add_switch(&log_frequency,
            NULL,
            "--log-frequency",
            "Frequency with which to log tree processing progress.");
    parser.add_switch(&quiet, "-q", "--quiet", "Suppress all informational/progress messages.");
    parser.parse(argc, argv);

    if (!quiet) {
        // std::string border(prog_id.size(), '=');
        std::string border(78, '=');
        colugo::console::out_ln(border);
        colugo::console::out_ln(prog_name, " - Phylogenetic STRUctural Distance EvaLuation (on Trees)");
        colugo::console::out_ln("[Revision: ", pstrudel::get_program_source_identification(), "]");
        colugo::console::out_ln("By Jeet Sukumaran");
        parser.write_description(std::cout);
        std::cout << "\n";
        colugo::console::out_ln(border);
    }

    if (!calculate_reference_distances && !calculate_pairwise_distances) {
        colugo::console::err_line("Need to specify at least one of '-r'/'--reference' or '-p'/'--pairwise' options");
        exit(EXIT_FAILURE);
    }

    // set up output filepaths
    if (COLUGO_FILESYS_PATH_SEPARATOR[0] == output_prefix[output_prefix.size()-1]) {
        output_prefix = default_output_filename_stem;
    }
    if (output_prefix[output_prefix.size() - 1] != '.') {
        output_prefix = output_prefix + '.';
    }
    std::map<std::string, std::string> output_filepaths;
    output_filepaths["log"] = output_prefix + "log";
    if (create_aggregated_comparison_trees_copy) {
        output_filepaths["comparison-trees"] = output_prefix + "comparison.trees";
    }
    if (calculate_reference_distances) {
        output_filepaths["reference-trees"] = output_prefix + "reference.trees";
        output_filepaths["reference-distances"] = output_prefix + "reference.distances.txt";
        output_filepaths["reference-distances-stacked"] = output_prefix + "reference.distances.stacked.txt";
    }
    if (calculate_pairwise_distances) {
        output_filepaths["pairwise-distances"] = output_prefix + "pairwise.distances.txt";
        output_filepaths["pairwise-distances-stacked"] = output_prefix + "pairwise.distances.stacked.txt";
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
            colugo::console::err_wrapped("\nRe-run the command using the '-x'/'--replace-existing-output' option to overwrite the files "
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

    // get and identify comparison trees
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

    // reference distances
    if (calculate_reference_distances) {
        // set up table
        platypus::DataTable results_table;
        results_table.add_key_column<unsigned long>("tree.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("source.filepath");
        }
        results_table.add_key_column<unsigned long>("num.tips");
        results_table.add_key_column<double>("tree.length");
        if (!reference_trees_filepath.empty()) {
            // user-supplied ref trees
            colugo::console::abort("User-specified reference trees not yet implemented");
        } else {
            // canonical ref trees
            results_table.add_data_column<double>("y.uw.max.unbalanced");
            results_table.add_data_column<double>("y.uw.max.unbalanced.scaled");
            results_table.add_data_column<double>("y.uw.max.balanced");
            results_table.add_data_column<double>("y.uw.max.balanced.scaled");
            if (calculate_symmetric_diff) {
                results_table.add_data_column<double>("urf.uw.max.unbalanced");
                results_table.add_data_column<double>("urf.uw.max.unbalanced.scaled");
                results_table.add_data_column<double>("urf.uw.max.balanced");
                results_table.add_data_column<double>("urf.uw.max.balanced.scaled");
            }
            TreePatternCollectionType tree_patterns;
            TreePatternMaxDistances tree_pattern_max_unweighted_pairwise_tip_profile_distances;
            TreePatternMaxDistances tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances;
            if (log_frequency == 0) {
                log_frequency = std::max(static_cast<unsigned long>(comparison_trees.size() / 10), 10UL);
            }

            logger.info("Begining calculating distances between comparison trees and canonical tree patterns");
            unsigned long comparison_tree_idx = 0;
            for (auto & comparison_tree : comparison_trees) {
                auto & results_table_row = results_table.add_row();

                auto comparison_tree_size = comparison_tree.get_num_tips();
                if (tree_patterns.find(comparison_tree_size) == tree_patterns.end()) {
                    logger.info("Building canonical ", comparison_tree_size, "-leaf reference trees");
                    build_canonical_tree_patterns(
                            comparison_tree_size,
                            tree_patterns,
                            tree_pattern_max_unweighted_pairwise_tip_profile_distances,
                            tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances,
                            true);
                }

                results_table_row.set("tree.idx", comparison_tree_idx + 1);
                if (add_tree_source_key) {
                    results_table_row.set("source.filepath", comparison_tree_sources.at(comparison_tree_idx));
                }
                results_table_row.set("num.tips", comparison_tree_size);
                results_table_row.set("tree.length", comparison_tree.get_total_tree_length());

                double d = 0.0;

                d = comparison_tree.get_unweighted_pairwise_tip_profile_distance(tree_patterns[comparison_tree_size]["max.unbalanced"]);
                results_table_row.set("y.uw.max.unbalanced", d);
                results_table_row.set("y.uw.max.unbalanced.scaled", d/tree_pattern_max_unweighted_pairwise_tip_profile_distances[comparison_tree_size]["max.unbalanced"]);

                d = comparison_tree.get_unweighted_pairwise_tip_profile_distance(tree_patterns[comparison_tree_size]["max.balanced"]);
                results_table_row.set("y.uw.max.balanced", d);
                results_table_row.set("y.uw.max.balanced.scaled", d/tree_pattern_max_unweighted_pairwise_tip_profile_distances[comparison_tree_size]["max.balanced"]);

                if (calculate_symmetric_diff) {
                    d = comparison_tree.get_unlabeled_symmetric_difference(tree_patterns[comparison_tree_size]["max.unbalanced"]);
                    results_table_row.set("urf.uw.max.unbalanced", d);
                    results_table_row.set("urf.uw.max.unbalanced.scaled", d/tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances[comparison_tree_size]["max.unbalanced"]);
                    d = comparison_tree.get_unlabeled_symmetric_difference(tree_patterns[comparison_tree_size]["max.balanced"]);
                    results_table_row.set("urf.uw.max.balanced", d);
                    results_table_row.set("urf.uw.max.balanced.scaled", d/tree_pattern_max_unweighted_unlabeled_symmetric_difference_distances[comparison_tree_size]["max.balanced"]);
                }

                comparison_tree_idx += 1;
            } // tree comparison

            // output canonical reference trees
            {
                auto & ref_trees_output_fpath = output_filepaths["reference-trees"];
                platypus::NewickWriter<pstrudel::DistanceTree> ref_tree_writer;
                ref_tree_writer.set_suppress_edge_lengths(true);
                platypus::bind_standard_interface(ref_tree_writer);
                std::ofstream ref_trees_out(ref_trees_output_fpath);
                for (auto & tp_by_size : tree_patterns) {
                    for (auto & tp_by_type : tp_by_size.second) {
                        ref_tree_writer.write(ref_trees_out, tp_by_type.second);
                        ref_trees_out << std::endl;
                    }
                }
            }

            logger.info("Completed calculating distances between comparison trees and canoncial tree patterns");

            // output primary results
            {
                auto & out_fpath = output_filepaths["reference-distances"];
                std::ofstream out(out_fpath);
                results_table.write(out);
            }

            // output stacked results
            {
                auto & out_fpath = output_filepaths["reference-distances-stacked"];
                std::ofstream out(out_fpath);
                results_table.write_stacked(out,
                        "dist.type", "dist");
            }

        } // canonical ref trees
    } // reference distances

    if (calculate_pairwise_distances) {
        // set up table
        platypus::DataTable results_table;
        results_table.add_key_column<unsigned long>("tree1.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree1.source.filepath");
        }
        results_table.add_key_column<unsigned long>("tree1.num.tips");
        results_table.add_key_column<double>("tree1.length");
        results_table.add_key_column<unsigned long>("tree2.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree2.source.filepath");
        }
        results_table.add_key_column<unsigned long>("tree2.num.tips");
        results_table.add_key_column<double>("tree2.length");

        results_table.add_data_column<double>("y.uw");
        results_table.add_data_column<double>("y.uw.scaled");
        results_table.add_data_column<double>("y.wt");
        results_table.add_data_column<double>("y.wt.scaled");
        if (calculate_symmetric_diff) {
            results_table.add_data_column<double>("urf.uw");
            results_table.add_data_column<double>("urf.uw.scaled");
        }
        logger.info("Begining calculating distances between all distinct pairs of trees");
        unsigned long num_trees = comparison_trees.size();
        unsigned long total_comparisons = num_trees * (num_trees - 1) / 2;
        unsigned long comparison_count = 0;
        for (unsigned long tree_idx1 = 0; tree_idx1 < num_trees - 1; ++tree_idx1) {
            auto & tree1 = comparison_trees[tree_idx1];
            for (unsigned long tree_idx2 = tree_idx1+1; tree_idx2 < num_trees; ++tree_idx2) {
                ++comparison_count;
                if (log_frequency > 0 && (comparison_count % log_frequency == 0)) {
                    logger.info("Comparison ", comparison_count, " of ", total_comparisons, ": Tree ", tree_idx1 + 1, " vs. tree ", tree_idx2 + 1);
                }
                auto & tree2 = comparison_trees[tree_idx2];
                auto & results_table_row = results_table.add_row();

                // tree 1 key
                results_table_row.set("tree1.idx", tree_idx1 + 1);
                if (add_tree_source_key) {
                    results_table_row.set("tree1.source.filepath", comparison_tree_sources.at(tree_idx1));
                }
                results_table_row.set("tree1.num.tips", tree1.get_num_tips());
                results_table_row.set("tree1.length", tree1.get_total_tree_length());

                // tree 2 key
                results_table_row.set("tree2.idx", tree_idx2 + 1);
                if (add_tree_source_key) {
                    results_table_row.set("tree2.source.filepath", comparison_tree_sources.at(tree_idx2));
                }
                results_table_row.set("tree2.num.tips", tree2.get_num_tips());
                results_table_row.set("tree2.length", tree2.get_total_tree_length());

                // data
                results_table_row.set("y.uw",
                        tree1.get_unweighted_pairwise_tip_profile_distance(tree2));
                // results_table_row.set("y.uw.scaled",
                //         tree1.get_unweighted_pairwise_tip_profile_distance(tree2));
                results_table_row.set("y.wt",
                        tree1.get_weighted_pairwise_tip_profile_distance(tree2));
                // results_table_row.set("y.wt.scaled",
                //         tree1.get_unweighted_pairwise_tip_profile_distance(tree2));
                if (calculate_symmetric_diff) {
                    results_table_row.set("urf.uw",
                            tree1.get_unlabeled_symmetric_difference(tree2));
                    // results_table_row.set("urf.uw.scaled",
                    //         tree1.get_unweighted_pairwise_tip_profile_distance(tree2));
                }
            }
        } // pairwise tree comparison

        // output primary results
        {
            auto & out_fpath = output_filepaths["pairwise-distances"];
            std::ofstream out(out_fpath);
            results_table.write(out);
        }

        // output stacked results
        {
            auto & out_fpath = output_filepaths["pairwise-distances-stacked"];
            std::ofstream out(out_fpath);
            results_table.write_stacked(out,
                    "dist.type", "dist");
        }
    } // pairwise distances
}


