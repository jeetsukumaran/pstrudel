#include <algorithm>
#include <cassert>
#include <colugo/filesys.hpp>
#include <colugo/cmdopt.hpp>
#include <colugo/logger.hpp>
#include <platypus/model/datatable.hpp>
#include <platypus/serialize/newick.hpp>
#include <platypus/parse/nclreader.hpp>
#include <platypus/model/treepattern.hpp>
#include <platypus/model/standardinterface.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pstrudel.hpp"
#include "distancetree.hpp"

//////////////////////////////////////////////////////////////////////////////
// WorkingTree

class WorkingTree : public pstrudel::DistanceTree {

    public:
        WorkingTree(const std::string & filepath, unsigned long file_index=0)
            : filepath_(filepath), file_index_(file_index), file_tree_index_(0) { }

        const std::string & get_filepath() const {
            return this->filepath_;
        }
        unsigned long get_file_index() const {
            return this->file_index_;
        }
        void set_file_index(unsigned long file_index) {
            this->file_index_ = file_index;
        }
        unsigned long get_file_tree_index() const {
            return this->file_tree_index_;
        }
        void set_file_tree_index(unsigned long file_tree_index) {
            this->file_tree_index_ = file_tree_index;
        }

    public:

    private:
        const std::string      filepath_;
        unsigned long          file_index_;
        unsigned long          file_tree_index_;

}; // WorkingTree

template <class TreeT>
void postprocess_working_tree(TreeT & tree, unsigned long idx, unsigned long ntips, unsigned long nints, double tree_length) {
    tree.set_file_tree_index(idx);
    tree.set_num_tips(ntips);
    tree.set_total_tree_length(tree_length);
}

template <class TreeT>
int get_trees(std::vector<TreeT>& trees,
        std::istream & src,
        const std::string & source_name,
        const std::string& format,
        unsigned long file_index) {
    auto reader = platypus::NclTreeReader<TreeT>();
    platypus::bind_standard_interface(reader);
    reader.set_tree_postprocess_fn(postprocess_working_tree<TreeT>);
    std::function<TreeT& ()> get_new_tree_reference = [&trees, &source_name, &file_index] () -> TreeT& { trees.emplace_back(source_name, file_index); return trees.back(); };
    return reader.read(src, get_new_tree_reference, format);
}

template <class TreeT>
    int get_trees(
            std::vector<TreeT>& trees,
            const std::vector<std::string> filepaths,
            const std::string& format,
            int type_of_file_info_to_store,
            colugo::Logger & logger) {
        auto reader = platypus::NclTreeReader<TreeT>();
        platypus::bind_standard_interface(reader);
        reader.set_tree_postprocess_fn(postprocess_working_tree<TreeT>);
        unsigned long file_idx = 0;
        unsigned long num_trees_read = 0;
        for (auto & filepath : filepaths) {
            auto fullpath = colugo::filesys::absolute_path(filepath);
            std::string stored_filepath;
            if (type_of_file_info_to_store <= 1) {
                stored_filepath = std::to_string(file_idx + 1);
            } else if (type_of_file_info_to_store <= 2) {
                stored_filepath = colugo::filesys::get_path_leaf(fullpath);
            } else if (type_of_file_info_to_store <= 3) {
                stored_filepath = filepath;
            } else if (type_of_file_info_to_store <= 4) {
                stored_filepath = fullpath;
            }
            logger.info("Reading trees for comparison from file ", file_idx+1, " of ", filepaths.size(), ": '", fullpath, "'");
            std::function<TreeT& ()> get_new_tree_reference = [&trees, &stored_filepath, &file_idx] () -> TreeT& { trees.emplace_back(stored_filepath, file_idx); return trees.back(); };
            std::ifstream src(fullpath);
            if (!src.good()) {
                logger.abort("Failed to open file for input: '", fullpath, "'");
            }
            num_trees_read += reader.read(src, get_new_tree_reference, format);
            // logger.info(num_trees_read, " trees read from '", fullpath,
            //         "' (total number of trees in comparison set is now: ", comparison_trees.size(), ")");
            ++file_idx;
        }
        return num_trees_read;
    }

    //////////////////////////////////////////////////////////////////////////////
    // TreePatternManager

    class CanonicalTreePatterns {
        public:
            typedef pstrudel::DistanceTree                         tree_type;
            typedef typename pstrudel::DistanceTree::value_type    value_type;
            CanonicalTreePatterns()
                    : num_tips_(0) {
            }
            std::vector<value_type> generate_leaves() {
                std::vector<value_type> leaves;
                leaves.reserve(this->num_tips_);
                for (unsigned long i = 0; i < this->num_tips_; ++i) {
                    leaves.emplace_back();
                    auto & leaf = leaves.back();
                    leaf.set_label("t" + std::to_string(i+1));
                }
                return leaves;
            }
            void generate(unsigned long num_tips) {
                // setup
                this->num_tips_ = num_tips;
                auto leaves = this->generate_leaves();

                // unbalanced tree, mean coalescent
                platypus::build_maximally_unbalanced_tree(this->tree_patterns_[".unbalanced.mean.coal"],
                        leaves.begin(), leaves.end());
                this->tree_patterns_[".unbalanced.mean.coal"].add_edge_lengths(0); // mean

                // unbalanced tree, converse coalescent
                platypus::build_maximally_unbalanced_tree(this->tree_patterns_[".unbalanced.converse.coal"],
                        leaves.begin(), leaves.end());
                this->tree_patterns_[".unbalanced.converse.coal"].add_edge_lengths(1); // converse

                // unbalanced tree, uniform coalescent
                platypus::build_maximally_unbalanced_tree(this->tree_patterns_[".unbalanced.uniform.coal"],
                        leaves.begin(), leaves.end());
                this->tree_patterns_[".unbalanced.uniform.coal"].add_edge_lengths(2); // uniform

                // balanced tree
                platypus::build_maximally_balanced_tree(this->tree_patterns_[".balanced"],
                        leaves.begin(), leaves.end());

                // calculate inter-reference tree distances
                this->tree_pattern_cross_distances_.add_key_column<std::string>("pattern");
                for (auto & tree_pattern_name : CanonicalTreePatterns::tree_pattern_names_) {
                this->tree_pattern_cross_distances_.add_data_column<double>("pwtd.uw" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("pwtd.wt" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("pwtd.swt" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("ltt" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("lst" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("lst.s" + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<unsigned long>("rfdu.uw" + tree_pattern_name);
            }
            for (auto & tree_pattern_name1 : CanonicalTreePatterns::tree_pattern_names_) {
                auto & row = this->tree_pattern_cross_distances_.add_row();
                row.set("pattern", tree_pattern_name1);
                auto & tree1 = this->tree_patterns_[tree_pattern_name1];
                for (auto & tree_pattern_name2 : CanonicalTreePatterns::tree_pattern_names_) {
                    auto & tree2 = this->tree_patterns_[tree_pattern_name2];
                    row.set("pwtd.uw" + tree_pattern_name2, tree2.get_unweighted_pairwise_tip_profile_distance(tree1));
                    row.set("pwtd.wt" + tree_pattern_name2, tree2.get_scaled_weighted_pairwise_tip_profile_distance(tree1));
                    row.set("pwtd.swt" + tree_pattern_name2, tree2.get_weighted_pairwise_tip_profile_distance(tree1));
                    row.set("ltt" + tree_pattern_name2, tree2.get_lineage_accumulation_profile_distance(tree1));
                    row.set("lst" + tree_pattern_name2, tree2.get_lineage_splitting_time_profile_distance(tree1));
                    row.set("lst.s" + tree_pattern_name2, tree2.get_scaled_lineage_splitting_time_profile_distance(tree1));
                    row.set("rfdu.uw" + tree_pattern_name2, tree2.get_unlabeled_symmetric_difference(tree1));
                }
            }
        }
        const std::vector<std::string> & get_tree_pattern_names() const {
            return CanonicalTreePatterns::tree_pattern_names_;
        }
        template <class T, class R>
        void score(
                const std::string tree_pattern_name,
                T & other_tree,
                R & row,
                bool scale_by_tree_length,
                bool calculate_symmetric_diff) {
            auto tpi = this->tree_patterns_.find(tree_pattern_name);
            if (tpi == this->tree_patterns_.end()) {
                colugo::console::abort("Invalid tree pattern name: '", tree_pattern_name, "'");
            }
            auto & tree_pattern = tpi->second;
            row.set("target.tree.pattern", tree_pattern_name);
            tree_pattern.tabulate_distances(
                    "",
                    other_tree,
                    row,
                    scale_by_tree_length,
                    calculate_symmetric_diff);
        }
        template <class W>
        void write_trees(W & writer, std::ostream & out) {
            for (auto & tree_pattern_name : CanonicalTreePatterns::tree_pattern_names_) {
                writer.write(out, this->tree_patterns_[tree_pattern_name]);
                out << "\n";
            }
        }
        void write_cross_distance_table(std::ostream & out) {
            this->tree_pattern_cross_distances_.write(out);
        }
    public:
        static void add_canonical_columns(platypus::DataTable & table,
                platypus::stream::OutputStreamFormatters & col_formatters,
                bool calculate_symmetric_diff) {
            table.add_key_column<std::string>("target.tree.pattern");
            pstrudel::DistanceTree::add_results_data_columns(
                    "",
                    table,
                    col_formatters,
                    calculate_symmetric_diff);
        }

    private:
        unsigned long                                    num_tips_;
        std::map<std::string, pstrudel::DistanceTree>    tree_patterns_;
        platypus::DataTable                              tree_pattern_cross_distances_;
        // std::map<std::string, double>                    max_unweighted_pairwise_tip_profile_distance_;
        // std::map<std::string, double>                    max_weighted_pairwise_tip_profile_distance_;
        // std::map<std::string, double>                    max_unlabeled_symmetric_difference_distance_;
        static const std::vector<std::string>            tree_pattern_names_;

}; // CanonicalTreePatterns
const std::vector<std::string> CanonicalTreePatterns::tree_pattern_names_{
        ".unbalanced.mean.coal",
        ".unbalanced.converse.coal",
        ".unbalanced.uniform.coal",
        ".balanced",
}; // static cons ttree_pattern_names_

//////////////////////////////////////////////////////////////////////////////
// main

int main(int argc, const char * argv[]) {
    std::string      prog_name                               = "pstrudel-trees";
    std::string      prog_id                                 = pstrudel::get_program_identification(prog_name).c_str();
    std::string      format                                  = "nexus";
    unsigned         long num_interpolated_points            = 0;
    bool             calculate_canonical_distances           = false;
    bool             calculate_pairwise_distances            = false;
    std::string      target_trees_filepath                   = "";
    bool             no_scale_by_tree_length                 = false;
    bool             calculate_symmetric_diff                = false;
    std::string      default_output_filename_stem            = "pstrudel-trees";
    std::string      output_prefix                           = default_output_filename_stem;
    bool             create_aggregated_comparison_trees_copy = false;
    bool             suppress_header_row                     = false;
    int              add_tree_source_key                     = 0;
    std::string      analysis_label;
    bool             replace_existing_output_files           = false;
    unsigned long    log_frequency                           = 0;
    bool             quiet                                   = false;

    colugo::OptionParser parser = colugo::OptionParser(
            prog_id.c_str(),
            "Calculate structural distances between unlabeled phylogenetic trees of various sizes.",
            "%prog [options] [TREE-FILE [TREE-FILE [...]]]");

    parser.add_switch(&calculate_pairwise_distances,
            "-p",
            "--pairwise",
            "Calculate distances from every tree to every other tree in input set.",
            nullptr,
            "Comparison Options");
    parser.add_switch(&calculate_canonical_distances,
            "-c",
            "--canonical",
            "Calculate distances from every tree in input set to canonical tree patterns.",
            nullptr,
            "Comparison Options");
    parser.add_option<std::string>(&target_trees_filepath, "-t", "--target-trees",
            "Calculate distances from every tree in input set to every tree in TARGET-TREE-FILEPATH.",
            "TARGET-TREE-FILEPATH",
            "Comparison Options");

    parser.add_switch(&no_scale_by_tree_length,
            nullptr,
            "--no-scale-by-tree-length",
            "By default, all edge weights will be scaled or normalized to tree length (sum of edge lengths on a tree) before distances are calculated;"
            " set this option to use unscaled (raw) edge lengths.",
            nullptr,
            "Distance Options");
    parser.add_switch(&calculate_symmetric_diff,
            "-d",
            "--calculate-symmetric-difference",
            "In addition to `y` or profile distances, calculate the Robinson-Foulds or symmetric difference.",
            nullptr,
            "Distance Options");

    parser.add_option<std::string>(&format, "-f", "--format",
            "Format for tree sources, one of: 'nexus' [default] or 'newick'.",
            "FORMAT",
            "Source Options");
    // parser.add_option<unsigned long>(&num_interpolated_points, "-n", "--profile-size",
    //         "Number of interpolated points in profile metric; if specified."
    //         " This must be equal or greater to than the largest input data size."
    //         " If not specified, will default to the largest input data size.",
    //         "SIZE");
    parser.add_option<std::string>(&output_prefix, "-o", "--output-prefix",
            "Prefix (directory path and filename stem) for output file(s); if"
            " not specified, will default to '%default'.",
            "PATH/TO/OUTPUT",
            "Output Options");
    parser.add_switch(&replace_existing_output_files,
            "-x",
            "--replace-existing-output",
            "Replace (overwrite) existing output files.",
            nullptr,
            "Output Options"
            );
    parser.add_switch(&add_tree_source_key, nullptr, "--add-tree-source-key",
            "Add a column in the results identifying the filename of the source of the tree(s) being compared:"
            " 1: value is index of file; 2: value is basename of file; 3: value is relative path of file;"
            " 4: value is absolute path of file.",
            nullptr, "Output Options");
    parser.add_option<std::string>(&analysis_label, "-l", "--add-label",
            "Add column to results with LABEL as values for all rows.",
            "LABEL",
            "Output Options");
    parser.add_switch(&suppress_header_row, nullptr, "--suppress-header-row",
            "Do not write column/field name row in results.",
            nullptr, "Output Options");
    parser.add_switch(&create_aggregated_comparison_trees_copy, nullptr, "--save-comparison-trees",
            "Save a copy of the comparison trees in a file (aggregating them from across multiple files if multiple source files specified).",
            nullptr, "Output Options");
    parser.add_switch(&log_frequency,
            nullptr,
            "--log-frequency",
            "Frequency with which to log tree processing progress.",
            nullptr, "Run Options");
    parser.add_switch(&quiet, "-q", "--quiet", "Suppress all informational/progress messages.",
            nullptr, "Run Options");
    parser.parse(argc, argv);
    bool scale_by_tree_length = !no_scale_by_tree_length;

    if (!quiet) {
        // std::string border(prog_id.size(), '=');
        std::string border(78, '=');
        std::string subborder(78, '-');
        colugo::console::out_ln(border);
        colugo::console::out_ln(prog_name, " - Phylogenetic STRUctural Distance EvaLuation (on Trees)");
        colugo::console::out_ln(subborder);
        colugo::console::out_ln("[Revision: ", pstrudel::get_program_source_identification(), "]");
        colugo::console::out_ln("By Jeet Sukumaran");
        parser.write_description(std::cout);
        std::cout << "\n";
        colugo::console::out_ln(border);
    }

    if (!calculate_canonical_distances && target_trees_filepath.empty() && !calculate_pairwise_distances) {
        colugo::console::err_line("Need to specify at least one of '-p'/'--pairwise', '-c'/'--canonical', or '-t'/'--target-trees' options");
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
    if (calculate_canonical_distances) {
        output_filepaths["canonical-distances"] = output_prefix + "canonical.distances.txt";
        output_filepaths["canonical-distances-stacked"] = output_prefix + "canonical.distances.stacked.txt";
    }
    if (!target_trees_filepath.empty()) {
        output_filepaths["target-distances"] = output_prefix + "target.distances.txt";
        output_filepaths["target-distances-stacked"] = output_prefix + "target.distances.stacked.txt";
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
    if (!logfile.good()) {
        colugo::console::abort("Unable to open log file for output: '", output_filepaths["log"], "'");
    }
    logger.add_channel(logfile,
            colugo::Logger::LoggingLevel::INFO,
            true,
            colugo::Logger::LoggingLevel::NOTSET);

    // log program identification
    logger.info(prog_id, " is running");

    // get and identify comparison trees
    std::vector<std::string> args = parser.get_args();
    std::vector<WorkingTree> comparison_trees;
    if (args.size() == 0) {
        logger.info("(reading trees for comparison from standard input)");
        get_trees(comparison_trees, std::cin, "<stdin>", format, 0);
        logger.info(comparison_trees.size(), " trees read from standard input");
    } else {
        get_trees(comparison_trees, args, format, add_tree_source_key, logger);
        logger.info(comparison_trees.size(), " trees read from ", args.size(), " file(s)");
    }

    platypus::stream::OutputStreamFormatters col_formatting{std::fixed, std::setprecision(16)};

    // target distances
    if (!target_trees_filepath.empty()) {
        std::vector<WorkingTree> target_trees;
        auto target_trees_full_filepath = colugo::filesys::absolute_path(target_trees_filepath);
        logger.info("Reading target trees from '", target_trees_full_filepath, "'");
        std::ifstream src(target_trees_full_filepath);
        get_trees(target_trees, src, target_trees_full_filepath, format, 0);
        logger.info(target_trees.size(), " trees read from '", target_trees_full_filepath, "'");
        platypus::DataTable results_table;
        if (!analysis_label.empty()) {
            results_table.add_key_column<std::string>("analysis");
        }
        results_table.add_key_column<unsigned long>("tree.i.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree.i.source.file");
            results_table.add_key_column<unsigned long>("tree.i.source.tree");
        }

        results_table.add_key_column<unsigned long>("target.tree.idx");
        pstrudel::DistanceTree::add_results_data_columns(
                "",
                results_table,
                col_formatting,
                calculate_symmetric_diff);

        // // each target tree gets its own column
        // unsigned long target_tree_idx = 0;
        // for (auto & ttree : target_trees) {
        //     std::string target_tree_label = ".t" + std::to_string(target_tree_idx+1);
        //     pstrudel::DistanceTree::add_results_data_columns(
        //             target_tree_label,
        //             results_table,
        //             col_formatting,
        //             calculate_symmetric_diff);
        //     ++target_tree_idx;
        // }

        logger.info("Beginning calculating distances between comparison trees and target tree(s)");
        if (scale_by_tree_length) {
            logger.info("Edge lengths will be scaled to tree length");
        } else {
            logger.info("Edge lengths will NOT be scaled to tree length");
        }
        unsigned long num_target_trees = target_trees.size();
        unsigned long num_comparison_trees = comparison_trees.size();
        unsigned long total_comparisons = num_comparison_trees * (num_comparison_trees - 1) / 2;
        unsigned long comparison_count = 0;
        unsigned long comparison_tree_idx = 0;
        unsigned long target_tree_idx = 0;
        for (auto & ctree : comparison_trees) {
            target_tree_idx = 0;
            for (auto & ttree : target_trees) {
                auto & results_table_row = results_table.add_row();
                if (log_frequency > 0 && (comparison_count % log_frequency == 0)) {
                    logger.info("Comparison ", comparison_count, " of ", total_comparisons, ": Target tree ", target_tree_idx + 1, " vs. tree ", comparison_tree_idx + 1);
                }
                if (!analysis_label.empty()) {
                    results_table_row.set("analysis", analysis_label);
                }
                results_table_row.set("tree.i.idx", comparison_tree_idx);
                if (add_tree_source_key) {
                    results_table_row.set("tree.i.source.file", ctree.get_filepath());
                    results_table_row.set("tree.i.source.tree", ctree.get_file_tree_index()+1);
                }
                results_table_row.set("target.tree.idx", target_tree_idx+1);
                ctree.tabulate_distances(
                        "",
                        ttree,
                        results_table_row,
                        scale_by_tree_length,
                        calculate_symmetric_diff);
                ++target_tree_idx;
                ++comparison_count;
            } // for each target tree
            ++comparison_tree_idx;
        } // for each tree in comparison set

        // output_filepaths["target-distances"] = output_prefix + "target.distances.txt";
        // output_filepaths["target-distances-stacked"] = output_prefix + "target.distances.stacked.txt";
        // output primary results
        {
            auto & out_fpath = output_filepaths["target-distances"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write(out, "\t", !suppress_header_row);
        }

        // output stacked results
        {
            auto & out_fpath = output_filepaths["target-distances-stacked"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write_stacked(out,
                    "dist.type", "dist", "\t", !suppress_header_row);
        }

    } // target distances

    // canonical distances
    if (calculate_canonical_distances) {
        // set up table
        platypus::DataTable results_table;
        if (!analysis_label.empty()) {
            results_table.add_key_column<std::string>("analysis");
        }
        results_table.add_key_column<unsigned long>("tree.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("source.file");
            results_table.add_key_column<std::string>("source.tree");
        }
        results_table.add_key_column<unsigned long>("num.tips");
        results_table.add_key_column<double>("tree.length", col_formatting);
        if (log_frequency == 0) {
            log_frequency = std::max(static_cast<unsigned long>(comparison_trees.size() / 10), 10UL);
        }
        logger.info("Beginning calculating distances between comparison trees and canonical tree patterns");
        if (scale_by_tree_length) {
            logger.info("Edge lengths will be scaled to tree length");
        } else {
            logger.info("Edge lengths will NOT be scaled to tree length");
        }
        std::map<unsigned long, CanonicalTreePatterns> tree_patterns;
        CanonicalTreePatterns::add_canonical_columns(results_table, col_formatting, calculate_symmetric_diff);
        unsigned long comparison_tree_idx = 0;
        for (auto & comparison_tree : comparison_trees) { // tree comparison
            auto comparison_tree_size = comparison_tree.get_num_tips();
            if (log_frequency > 0 && (comparison_tree_idx % log_frequency == 0)) {
                logger.info("Calculating canonical distances for tree ", comparison_tree_idx + 1, " of ", comparison_tree_size);
            }
            if (tree_patterns.find(comparison_tree_size) == tree_patterns.end()) {
                logger.info("Building canonical ", comparison_tree_size, "-leaf reference trees");
                tree_patterns[comparison_tree_size].generate(comparison_tree_size);
            }
            for (auto & tree_pattern_name : tree_patterns[comparison_tree_size].get_tree_pattern_names()) {
                auto & results_table_row = results_table.add_row();
                if (!analysis_label.empty()) {
                    results_table_row.set("analysis", analysis_label);
                }
                results_table_row.set("tree.idx", comparison_tree_idx + 1);
                if (add_tree_source_key) {
                    results_table_row.set("source.file", comparison_tree.get_filepath());
                    results_table_row.set("source.tree", comparison_tree.get_file_tree_index()+1);
                }
                results_table_row.set("num.tips", comparison_tree_size);
                results_table_row.set("tree.length", comparison_tree.get_total_tree_length());
                tree_patterns[comparison_tree_size].score(tree_pattern_name,
                        comparison_tree,
                        results_table_row,
                        scale_by_tree_length,
                        calculate_symmetric_diff);
            }
            comparison_tree_idx += 1;
        } // tree comparison
        logger.info("Completed calculating distances between comparison trees and canoncial tree patterns");

        // output canonical info
        {
            platypus::NewickWriter<pstrudel::DistanceTree> ref_tree_writer;
            ref_tree_writer.set_suppress_edge_lengths(false);
            platypus::bind_standard_interface(ref_tree_writer);
            for (auto & tpi : tree_patterns) {
                std::string prefix = output_prefix + "canonical.n" + std::to_string(tpi.first) + ".";
                std::ofstream trees_out(prefix + "trees");
                tpi.second.write_trees(ref_tree_writer, trees_out);
                std::ofstream table_out(prefix + "distances.txt");
                tpi.second.write_cross_distance_table(table_out);
            }
        }

        // output primary results
        {
            auto & out_fpath = output_filepaths["canonical-distances"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write(out, "\t", !suppress_header_row);
        }

        // output stacked results
        {
            auto & out_fpath = output_filepaths["canonical-distances-stacked"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write_stacked(out,
                    "dist.type", "dist", "\t", !suppress_header_row);
        }
    } // canonical ref trees

    // pairwise distances
    if (calculate_pairwise_distances) {
        // set up table
        platypus::DataTable results_table;
        if (!analysis_label.empty()) {
            results_table.add_key_column<std::string>("analysis");
        }
        results_table.add_key_column<unsigned long>("tree.i.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree.i.source.file");
            results_table.add_key_column<unsigned long>("tree.i.source.tree");
        }
        results_table.add_key_column<unsigned long>("tree.i.num.tips");
        results_table.add_key_column<double>("tree.i.length", col_formatting);
        results_table.add_key_column<unsigned long>("tree.j.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree.j.source.file");
            results_table.add_key_column<unsigned long>("tree.j.source.tree");
        }
        results_table.add_key_column<unsigned long>("tree.j.num.tips");
        results_table.add_key_column<double>("tree.j.length", col_formatting);

        results_table.add_key_column<std::string>("comp.type");
        results_table.add_key_column<std::string>("comp.class");

        pstrudel::DistanceTree::add_results_data_columns(
                "",
                results_table,
                col_formatting,
                calculate_symmetric_diff);

        pstrudel::DistanceTree::add_results_data_columns(
                ".norm",
                results_table,
                col_formatting,
                calculate_symmetric_diff);

        logger.info("Beginning calculating distances between all distinct pairs of trees");
        if (scale_by_tree_length) {
            logger.info("Edge lengths will be scaled to tree length");
        } else {
            logger.info("Edge lengths will NOT be scaled to tree length");
        }
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

                // analysis id
                if (!analysis_label.empty()) {
                    results_table_row.set("analysis", analysis_label);
                }

                // tree 1 key
                results_table_row.set("tree.i.idx", tree_idx1 + 1);
                if (add_tree_source_key) {
                    results_table_row.set("tree.i.source.file", tree1.get_filepath());
                    results_table_row.set("tree.i.source.tree", tree1.get_file_tree_index()+1);
                }
                results_table_row.set("tree.i.num.tips", tree1.get_num_tips());
                results_table_row.set("tree.i.length", tree1.get_total_tree_length());

                // tree 2 key
                results_table_row.set("tree.j.idx", tree_idx2 + 1);
                if (add_tree_source_key) {
                    results_table_row.set("tree.j.source.file", tree2.get_filepath());
                    results_table_row.set("tree.j.source.tree", tree2.get_file_tree_index()+1);
                }
                results_table_row.set("tree.j.num.tips", tree2.get_num_tips());
                results_table_row.set("tree.j.length", tree2.get_total_tree_length());

                // comparison key
                if (tree1.get_filepath() == tree2.get_filepath()) {
                    results_table_row.set("comp.type", "within");
                } else {
                    results_table_row.set("comp.type", "between");
                }
                results_table_row.set("comp.class", std::to_string(tree1.get_file_index()+1) + ":" + std::to_string(tree2.get_file_index()+1));

                // data
                tree1.tabulate_distances(
                        "",
                        tree2,
                        results_table_row,
                        scale_by_tree_length,
                        calculate_symmetric_diff);

                ++comparison_count;
            } // for each tree_idx2 pairwise tree comparison
        } // for each tree_idx1 pairwise tree comparison

        logger.info("Calculating normalized distances");
        auto max_y_ptd_uw = results_table.column("pwtd.uw").max<double>();
        auto max_y_ptd_wt = results_table.column("pwtd").max<double>();
        auto max_y_ltt    = results_table.column("ltt").max<double>();
        auto max_y_lst    = results_table.column("lst").max<double>();
        double max_rfd_uw = 1.0;
        if (calculate_symmetric_diff) {
            max_rfd_uw = results_table.column("rfdu").max<double>();
        }
        unsigned long row_idx = 0;
        unsigned long nrows = results_table.num_rows();
        for (auto & row : results_table) {
            if (log_frequency > 0 && (row_idx % log_frequency == 0)) {
                logger.info("Normalizing entry ", row_idx + 1, " of ", nrows);
            }
            row.set("pwtd.uw.norm", row.get<double>("pwtd.uw") / max_y_ptd_uw);
            row.set("pwtd.norm", row.get<double>("pwtd") / max_y_ptd_wt);
            row.set("ltt.norm", row.get<double>("ltt") / max_y_ltt);
            row.set("lst.norm", row.get<double>("lst") / max_y_lst);
            if (calculate_symmetric_diff) {
                row.set("rfdu.norm", row.get<double>("rfdu") / max_rfd_uw);
            }
            ++row_idx;
        }
        logger.info("Completed calculating pairwise distances between all distinct pairs of trees");

        // output primary results
        {
            auto & out_fpath = output_filepaths["pairwise-distances"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write(out, "\t", !suppress_header_row);
        }

        // output stacked results
        {
            auto & out_fpath = output_filepaths["pairwise-distances-stacked"];
            std::ofstream out(out_fpath);
            if (!out.good()) {
                logger.abort("Unable to open file for output: '", out_fpath, "'");
            }
            results_table.write_stacked(out,
                    "dist.type", "dist", "\t", !suppress_header_row);
        }
    } // pairwise distances
}


