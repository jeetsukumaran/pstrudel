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
        WorkingTree(const std::string & filepath)
            : filepath_(filepath), file_index_(0) { }

        const std::string & get_filepath() const {
            return this->filepath_;
        }
        unsigned long get_file_index() const {
            return this->file_index_;
        }
        void set_file_index(unsigned long file_index) {
            this->file_index_ = file_index;
        }

    private:
        const std::string      filepath_;
        unsigned long          file_index_;


}; // WorkingTree

template <class TreeT>
void postprocess_working_tree(TreeT & tree, unsigned long idx, unsigned long ntips, unsigned long nints, double tree_length) {
    tree.set_file_index(idx);
    tree.set_num_tips(ntips);
    tree.set_total_tree_length(tree_length);
}

template <class TreeT>
int get_trees(std::vector<TreeT>& trees,
        const std::string & filepath,
        const std::string& format) {
    auto reader = platypus::NclTreeReader<TreeT>();
    platypus::bind_standard_interface(reader);
    reader.set_tree_postprocess_fn(postprocess_working_tree<TreeT>);
    std::function<TreeT& ()> get_new_tree_reference = [&trees, &filepath] () -> TreeT& { trees.emplace_back(filepath); return trees.back(); };
    return reader.read(std::ifstream(filepath), get_new_tree_reference, format);
}

template <class TreeT>
int get_trees(std::vector<TreeT>& trees,
        std::istream & src,
        const std::string & source_name,
        const std::string& format) {
    auto reader = platypus::NclTreeReader<TreeT>();
    platypus::bind_standard_interface(reader);
    reader.set_tree_postprocess_fn(postprocess_working_tree<TreeT>);
    std::function<TreeT& ()> get_new_tree_reference = [&trees, &source_name] () -> TreeT& { trees.emplace_back(source_name); return trees.back(); };
    return reader.read(src, get_new_tree_reference, format);
}

//////////////////////////////////////////////////////////////////////////////
// TreePatternManager

class CanonicalTreePatterns {
    public:
        typedef pstrudel::DistanceTree                         tree_type;
        typedef typename pstrudel::DistanceTree::value_type    value_type;
        CanonicalTreePatterns()
                : num_tips_(0)
                  , is_labeled_(true)
                  , calculate_symmetric_diff_(true) {
        }
        std::vector<value_type> generate_leaves() {
            std::vector<value_type> leaves;
            leaves.reserve(this->num_tips_);
            for (unsigned long i = 0; i < this->num_tips_; ++i) {
                leaves.emplace_back();
                if (this->is_labeled_) {
                    auto & leaf = leaves.back();
                    leaf.set_label("t" + std::to_string(i+1));
                }
            }
            return leaves;
        }
        void generate(unsigned long num_tips, bool is_labled=true, bool calculate_symmetric_diff=true) {
            this->num_tips_ = num_tips;
            this->calculate_symmetric_diff_ = calculate_symmetric_diff;
            auto leaves = this->generate_leaves();
            auto & t1 = this->tree_patterns_[CanonicalTreePatterns::tree_pattern_names_[0]];
            platypus::build_maximally_unbalanced_tree(t1, leaves.begin(), leaves.end());
            auto & t2 = this->tree_patterns_[CanonicalTreePatterns::tree_pattern_names_[1]];
            platypus::build_maximally_balanced_tree(t2, leaves.begin(), leaves.end());

            this->tree_pattern_cross_distances_.add_key_column<std::string>("pattern");
            for (auto & tree_pattern_name : CanonicalTreePatterns::tree_pattern_names_) {
                this->tree_pattern_cross_distances_.add_data_column<double>("y.uw." + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("y.wt." + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<double>("y.swt." + tree_pattern_name);
                this->tree_pattern_cross_distances_.add_data_column<unsigned long>("urf.uw." + tree_pattern_name);
            }

            for (auto & tree_pattern_name1 : CanonicalTreePatterns::tree_pattern_names_) {
                auto & row = this->tree_pattern_cross_distances_.add_row();
                row.set("pattern", tree_pattern_name1);
                auto & tree1 = this->tree_patterns_[tree_pattern_name1];
                for (auto & tree_pattern_name2 : CanonicalTreePatterns::tree_pattern_names_) {
                    auto & tree2 = this->tree_patterns_[tree_pattern_name2];
                    row.set ("y.uw." + tree_pattern_name2, tree2.get_unweighted_pairwise_tip_profile_distance(tree1));
                    row.set ("y.wt." + tree_pattern_name2, tree2.get_weighted_pairwise_tip_profile_distance(tree1));
                    row.set("y.swt." + tree_pattern_name2, tree2.get_scaled_weighted_pairwise_tip_profile_distance(tree1));
                    row.set("urf.uw." + tree_pattern_name2, tree2.get_unlabeled_symmetric_difference(tree1));
                }
            }

            for (auto & tree_pattern_name : CanonicalTreePatterns::tree_pattern_names_) {
                this->max_unweighted_pairwise_tip_profile_distance_[tree_pattern_name] = this->tree_pattern_cross_distances_.column("y.uw." + tree_pattern_name).max<double>();
                this->max_weighted_pairwise_tip_profile_distance_[tree_pattern_name] = this->tree_pattern_cross_distances_.column("y.wt." + tree_pattern_name).max<double>();
                this->max_scaled_weighted_pairwise_tip_profile_distance_[tree_pattern_name] = this->tree_pattern_cross_distances_.column("y.swt." + tree_pattern_name).max<double>();
                this->max_unlabeled_symmetric_difference_distance_[tree_pattern_name] = this->tree_pattern_cross_distances_.column("urf.uw." + tree_pattern_name).max<double>();
            }

        }
        template <class T, class R>
        void score(T & other_tree, R & row) {
            double d = 0.0;
            for (auto & tpi : this->tree_patterns_) {
                auto & tree_pattern_name = tpi.first;
                auto & tree_pattern = tpi.second;
                d = tree_pattern.get_unweighted_pairwise_tip_profile_distance(other_tree);
                row.set("y.uw." + tree_pattern_name, d);
                row.set("y.uw." + tree_pattern_name + ".scaled", d/this->max_unweighted_pairwise_tip_profile_distance_[tree_pattern_name]);

                d = tree_pattern.get_weighted_pairwise_tip_profile_distance(other_tree);
                row.set("y.wt." + tree_pattern_name, d);
                row.set("y.wt." + tree_pattern_name +".scaled", d/this->max_weighted_pairwise_tip_profile_distance_[tree_pattern_name]);

                d = tree_pattern.get_scaled_weighted_pairwise_tip_profile_distance(other_tree);
                row.set("y.swt." + tree_pattern_name, d);
                row.set("y.swt." + tree_pattern_name + ".scaled", d/this->max_scaled_weighted_pairwise_tip_profile_distance_[tree_pattern_name]);

                if (this->calculate_symmetric_diff_) {
                    d = tree_pattern.get_unlabeled_symmetric_difference(other_tree);
                    row.set("urf.uw." + tree_pattern_name, d);
                    row.set("urf.uw." + tree_pattern_name + ".scaled", d/this->max_unlabeled_symmetric_difference_distance_[tree_pattern_name]);
                }
            }
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
        static void add_results_data_columns(platypus::DataTable & table,
                platypus::stream::OutputStreamFormatters & col_formatters,
                bool calculate_symmetric_diff) {
            for (auto & tree_pattern_name : CanonicalTreePatterns::tree_pattern_names_) {
                for (auto & y_distance_name : CanonicalTreePatterns::tree_pattern_y_distance_names_) {
                    table.add_data_column<double>(y_distance_name + "." + tree_pattern_name, col_formatters);
                    table.add_data_column<double>(y_distance_name + "." + tree_pattern_name + ".scaled", col_formatters);
                }
                if (calculate_symmetric_diff) {
                    table.add_data_column<unsigned long>("urf.uw." + tree_pattern_name, col_formatters);
                    table.add_data_column<double>("urf.uw." + tree_pattern_name + ".scaled", col_formatters);
                }
            }
        }

    private:
        unsigned long                                    num_tips_;
        bool                                             is_labeled_;
        bool                                             calculate_symmetric_diff_;
        std::map<std::string, pstrudel::DistanceTree>    tree_patterns_;
        platypus::DataTable                              tree_pattern_cross_distances_;
        std::map<std::string, double>                    max_unweighted_pairwise_tip_profile_distance_;
        std::map<std::string, double>                    max_weighted_pairwise_tip_profile_distance_;
        std::map<std::string, double>                    max_scaled_weighted_pairwise_tip_profile_distance_;
        std::map<std::string, double>                    max_unlabeled_symmetric_difference_distance_;
        static const std::vector<std::string>            tree_pattern_names_;
        static const std::vector<std::string>            tree_pattern_y_distance_names_;

}; // CanonicalTreePatterns
const std::vector<std::string> CanonicalTreePatterns::tree_pattern_names_{
        "max.unbalanced",
        "max.balanced",
}; // static cons ttree_pattern_names_
const std::vector<std::string> CanonicalTreePatterns::tree_pattern_y_distance_names_{
        "y.uw",
        "y.wt",
        "y.swt",
}; // static cons ttree_pattern_names_

//////////////////////////////////////////////////////////////////////////////
// main

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
    std::vector<WorkingTree> comparison_trees;
    if (args.size() == 0) {
        logger.info("(reading trees for comparison from standard input)");
        get_trees(comparison_trees, std::cin, "<stdin>", format);
        logger.info(comparison_trees.size(), " trees read from standard input");
    } else {
        unsigned file_idx = 1;
        unsigned long num_trees_read = 0;
        for (auto & arg : args) {
            auto fullpath = colugo::filesys::absolute_path(arg);
            logger.info("Reading trees for comparison from file ", file_idx, " of ", args.size(), ": '", fullpath, "'");
            num_trees_read = get_trees(comparison_trees, fullpath, format);
            logger.info(num_trees_read, " trees read from '", fullpath,
                    "' (total number of trees in comparison set is now: ", comparison_trees.size(), ")");
            ++file_idx;
        }
        logger.info(comparison_trees.size(), " trees read from ", args.size(), " file(s)");
    }

    platypus::stream::OutputStreamFormatters col_formatting{std::fixed, std::setprecision(16)};

    // reference distances
    if (calculate_reference_distances) {
        // set up table
        platypus::DataTable results_table;
        results_table.add_key_column<unsigned long>("tree.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("source.filepath");
            results_table.add_key_column<std::string>("source.file.index");
        }
        results_table.add_key_column<unsigned long>("num.tips");
        results_table.add_key_column<double>("tree.length", col_formatting);
        if (!reference_trees_filepath.empty()) {
            // user-supplied ref trees
            colugo::console::abort("User-specified reference trees not yet implemented");
        } else {
            // canonical ref trees
            if (log_frequency == 0) {
                log_frequency = std::max(static_cast<unsigned long>(comparison_trees.size() / 10), 10UL);
            }
            logger.info("Beginning calculating distances between comparison trees and canonical tree patterns");
            std::map<unsigned long, CanonicalTreePatterns> tree_patterns;
            CanonicalTreePatterns::add_results_data_columns(results_table, col_formatting, calculate_symmetric_diff);
            unsigned long comparison_tree_idx = 0;
            for (auto & comparison_tree : comparison_trees) {
                auto comparison_tree_size = comparison_tree.get_num_tips();
                auto & results_table_row = results_table.add_row();
                results_table_row.set("tree.idx", comparison_tree_idx + 1);
                if (add_tree_source_key) {
                    results_table_row.set("source.filepath", comparison_tree.get_filepath());
                    results_table_row.set("source.file.index", comparison_tree.get_file_index()+1);
                }
                results_table_row.set("num.tips", comparison_tree_size);
                results_table_row.set("tree.length", comparison_tree.get_total_tree_length());

                if (tree_patterns.find(comparison_tree_size) == tree_patterns.end()) {
                    logger.info("Building canonical ", comparison_tree_size, "-leaf reference trees");
                    tree_patterns[comparison_tree_size].generate(comparison_tree_size, true, calculate_symmetric_diff);
                }
                tree_patterns[comparison_tree_size].score(comparison_tree, results_table_row);

                comparison_tree_idx += 1;
            } // tree comparison
            logger.info("Completed calculating distances between comparison trees and canoncial tree patterns");

            // output canonical info
            {
                platypus::NewickWriter<pstrudel::DistanceTree> ref_tree_writer;
                ref_tree_writer.set_suppress_edge_lengths(true);
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
            results_table.add_key_column<unsigned long>("tree1.source.file.index");
        }
        results_table.add_key_column<unsigned long>("tree1.num.tips");
        results_table.add_key_column<double>("tree1.length", col_formatting);
        results_table.add_key_column<unsigned long>("tree2.idx");
        if (add_tree_source_key) {
            results_table.add_key_column<std::string>("tree2.source.filepath");
            results_table.add_key_column<unsigned long>("tree2.source.file.index");
        }
        results_table.add_key_column<unsigned long>("tree2.num.tips");
        results_table.add_key_column<double>("tree2.length", col_formatting);

        results_table.add_data_column<double>("y.uw", col_formatting);
        results_table.add_data_column<double>("y.uw.norm", col_formatting);
        results_table.add_data_column<double>("y.wt", col_formatting);
        results_table.add_data_column<double>("y.wt.norm", col_formatting);
        results_table.add_data_column<double>("urf.uw");
        results_table.add_data_column<double>("urf.uw.norm", col_formatting);
        if (!calculate_symmetric_diff) {
            results_table.column("urf.uw").set_hidden(true);
            results_table.column("urf.uw.norm").set_hidden(true);
        }
        logger.info("Beginning calculating distances between all distinct pairs of trees");
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
                    results_table_row.set("tree1.source.filepath", tree1.get_filepath());
                    results_table_row.set("tree1.source.file.index", tree1.get_file_index()+1);
                }
                results_table_row.set("tree1.num.tips", tree1.get_num_tips());
                results_table_row.set("tree1.length", tree1.get_total_tree_length());

                // tree 2 key
                results_table_row.set("tree2.idx", tree_idx2 + 1);
                if (add_tree_source_key) {
                    results_table_row.set("tree2.source.filepath", tree2.get_filepath());
                    results_table_row.set("tree2.source.file.index", tree2.get_file_index()+1);
                }
                results_table_row.set("tree2.num.tips", tree2.get_num_tips());
                results_table_row.set("tree2.length", tree2.get_total_tree_length());

                // data
                results_table_row.set("y.uw",
                        tree1.get_unweighted_pairwise_tip_profile_distance(tree2));
                results_table_row.set("y.wt",
                        tree1.get_weighted_pairwise_tip_profile_distance(tree2));
                results_table_row.set("urf.uw",
                        tree1.get_unlabeled_symmetric_difference(tree2));
            }
        } // pairwise tree comparison

        logger.info("Calculating normalized distances");
        auto y_uw_col = results_table.get_column<double>("y.uw");
        auto y_wt_col = results_table.get_column<double>("y.wt");
        auto urf_uw_col = results_table.get_column<double>("urf.uw");
        auto max_y_uw = *(std::max_element(y_uw_col.begin(), y_uw_col.end()));
        auto max_y_wt = *(std::max_element(y_wt_col.begin(), y_wt_col.end()));
        auto max_urf_uw = *(std::max_element(urf_uw_col.begin(), urf_uw_col.end()));
        for (auto & row : results_table) {
            row.set("y.uw.norm", row.get<double>("y.uw") / max_y_uw);
            row.set("y.wt.norm", row.get<double>("y.wt") / max_y_wt);
            row.set("urf.uw.norm", row.get<double>("urf.uw") / max_urf_uw);
        }
        logger.info("Completed calculating pairwise distances between all distinct pairs of trees");

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


