#include <algorithm>
#include <cassert>
#include <colugo/filesys.hpp>
#include <colugo/cmdopt.hpp>
#include <colugo/logger.hpp>
#include <platypus/parse/nclreader.hpp>
#include <platypus/parse/newick.hpp>
#include <platypus/model/standardinterface.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pstrudel.hpp"
#include "treeshape.hpp"

const unsigned long DEFAULT_PRECISION = 22;

//////////////////////////////////////////////////////////////////////////////
// WorkingTree

class WorkingTree : public pstrudel::TreeShape {

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

//////////////////////////////////////////////////////////////////////////////
// Utility functions

template <class TreeT>
int get_trees(std::vector<TreeT>& trees,
        std::istream & src,
        const std::string & source_name,
        const std::string& format,
        unsigned long file_index) {
    std::function<TreeT& ()> get_new_tree_reference =
        [&trees, &source_name, &file_index] () -> TreeT& {
            trees.emplace_back(source_name, file_index); return trees.back();
        };
    if (format == "newick") {
        auto reader = platypus::NewickReader<TreeT>();
        platypus::bind_standard_interface(reader);
        return reader.read(src, get_new_tree_reference);
    } else {
        auto reader = platypus::NclTreeReader<TreeT>();
        platypus::bind_standard_interface(reader);
        return reader.read(src, get_new_tree_reference, format);
    }
}

//////////////////////////////////////////////////////////////////////////////
// main

int main(int argc, const char * argv[]) {
    std::string      prog_name                               = "pstrudel-tree-profile-calc";
    std::string      prog_id                                 = pstrudel::get_program_identification(prog_name).c_str();
    std::string      format                                  = "nexus";
    std::string      reference_trees_filepath                = "";
    bool             no_scale_by_tree_length                 = false;
    std::string      default_output_filename_stem            = "pstrudel-tree-profiles";
    std::string      output_prefix                           = default_output_filename_stem;
    unsigned long    output_precision                        = DEFAULT_PRECISION;
    std::string      analysis_label                          ;
    unsigned int     add_tree_source_key                     = 0;
    bool             replace_existing_output_files           = false;
    bool             quiet                                   = false;

    colugo::OptionParser parser = colugo::OptionParser(
            prog_id.c_str(),
            "Calculate phylogenetic tree profiles.",
            "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    parser.add_option<std::string>(&reference_trees_filepath, "-r", "--reference-trees",
            "Calculate distances from every tree in input set to every tree in REFERENCE-TREE-FILEPATH.",
            "REFERENCE-TREE-FILEPATH",
            "Comparison Options");
    parser.add_switch(&no_scale_by_tree_length,
            nullptr,
            "--no-scale-by-tree-length",
            "By default, for all distance calculations, edge weights will be scaled or normalized to tree length (sum of edge lengths on a tree) before distances are calculated;"
            " set this option to use unscaled (raw) edge lengths.",
            nullptr,
            "Calculation Options");
    parser.add_option<std::string>(&format, "-f", "--format",
            "Format for tree sources, one of: 'nexus' [default] or 'newick'.",
            "FORMAT",
            "Source Options");
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
    parser.add_option<unsigned long>(&output_precision, nullptr, "--output-precision",
            "Number of digits of precision in output (default=%default).",
            "NUM-DIGITS", "Output Options");
    parser.add_option<std::string>(&analysis_label, "-l", "--add-label",
            "Add column to results with LABEL as values for all rows.",
            "LABEL",
            "Output Options");
    parser.add_option<unsigned int>(&add_tree_source_key, nullptr, "--add-tree-source-key",
            "Add a column in the results identifying the filename of the source of the tree(s) being compared:"
            " 1: value is index of file; 2: value is basename of file; 3: value is relative path of file;"
            " 4: value is absolute path of file.",
            nullptr, "Output Options");
    parser.add_switch(&quiet, "-q", "--quiet", "Suppress all informational/progress messages.",
            nullptr, "Run Options");
    parser.parse(argc, argv);
    bool scale_by_tree_length = !no_scale_by_tree_length;

    if (!quiet) {
        // std::string border(prog_id.size(), '=');
        std::string border(78, '=');
        std::string subborder(78, '-');
        colugo::console::out_ln(border);
        colugo::console::out_ln(prog_name, " - Phylogenetic STRUctural Distance EvaLuation - Tree Profile Calculator");
        colugo::console::out_ln(subborder);
        colugo::console::out_ln("[Revision: ", pstrudel::get_program_source_identification(), "]");
        colugo::console::out_ln("By Jeet Sukumaran");
        parser.write_description(std::cout);
        std::cout << "\n";
        colugo::console::out_ln(border);
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
    output_filepaths["dists"] = output_prefix + "dists";

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

    // idiot checking
    std::vector<std::string> args = parser.get_args();
    if (args.empty()) {
        logger.abort("Source of comparison tree(s) not specified");
    }

    // comparison to reference trees
    if (reference_trees_filepath.empty()) {
        logger.abort("Source of reference tree(s) not specified");
    } else {

        {
        logger.info("Reading reference trees");
        std::vector<WorkingTree> reference_trees;
        auto reference_trees_full_filepath = colugo::filesys::absolute_path(reference_trees_filepath);
        logger.info("Reading reference trees from '", reference_trees_full_filepath, "'");
        std::ifstream src(reference_trees_full_filepath);
        get_trees(
                reference_trees,
                src,
                reference_trees_full_filepath,
                format,
                0);
        logger.info("'", reference_trees_full_filepath, "': ", reference_trees.size());
        }

        {
        logger.info("Reading comparison trees");
        std::vector<WorkingTree> comparison_trees;
        unsigned long file_idx = 0;
        unsigned long trees_read = 0;
        for (auto & arg : args) {
            auto comparison_trees_full_filepath = colugo::filesys::absolute_path(arg);
            logger.info("Reading comparison trees from comparison source file ",
                    file_idx+1,
                    " of ",
                    args.size(),
                    ": '",
                    comparison_trees_full_filepath,
                    "'");
            std::ifstream src(comparison_trees_full_filepath);
            trees_read = get_trees(
                    comparison_trees,
                    src,
                    comparison_trees_full_filepath,
                    format,
                    file_idx+1);
            logger.info(
                    "'",
                    comparison_trees_full_filepath,
                    "': ",
                    trees_read,
                    " trees read (total of ",
                    comparison_trees.size(),
                    " trees now in comparison set)");
            ++file_idx;
        }
        }

    }

}


