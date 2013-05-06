#include <cassert>
#include <colugo/filesys.hpp>
#include <colugo/cmdopt.hpp>
#include <colugo/logger.hpp>
#include "dataio.hpp"
#include "split.hpp"
#include "pstrudel.hpp"
#include "distancetree.hpp"

int main(int argc, const char * argv[]) {
    std::string      prog_name                            = "pstrudel-trees";
    std::string      prog_id                              = pstrudel::get_program_identification(prog_name).c_str();
    std::string      format                               = "nexus";
    unsigned         long num_interpolated_points         = 0;
    bool             calculate_reference_distances        = false;
    bool             calculate_pairwise_distances         = false;
    bool             calculate_symmetric_diff             = false;
    std::string      reference_trees_filepath             = "";
    std::string      default_output_filename_stem         = "pstrudel-trees";
    std::string      output_prefix                        = default_output_filename_stem;
    bool             suppress_copying_of_comparison_trees = false;
    bool             suppress_header_row                  = false;
    bool             replace_existing_output_files        = false;
    bool             quiet                                = false;

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
    parser.add_switch(&suppress_copying_of_comparison_trees, NULL, "--suppress-comparison-tree-copy",
            "Do not save an aggregated copy of the comparison tree files.");
    parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row",
            "Do not write column/field name row in results.");
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
    if (!suppress_copying_of_comparison_trees) {
        output_filepaths["comparison-trees"] = output_prefix + "comparison.trees";
    }
    if (calculate_reference_distances) {
        output_filepaths["reference-trees"] = output_prefix + "reference.trees";
        output_filepaths["reference-distances"] = output_prefix + "distances.reference.txt";
        output_filepaths["reference-distances-stacked"] = output_prefix + "distances.reference.stacked.txt";
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
}


