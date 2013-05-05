
int main(int argc, const char * argv[]) {
    // std::string prog_id = pstrudel::get_program_identification("PSTRUDEL-TREES-REFDIST").c_str();
    // std::string format = "nexus";
    // // bool force_rooted = false;
    // // bool force_unrooted = false;
    // unsigned long num_interpolated_points = 0;
    // bool exclude_distance_between_reference_trees = false;
    // bool calculate_other_distance_metrics = false;
    // std::string output_prefix = "";
    // bool suppress_header_row = false;
    // bool quiet = false;

    // colugo::OptionParser parser = colugo::OptionParser(
    //         prog_id.c_str(),
    //         "Calculate and report distances between each tree and reference trees (maximally-balanced and maximally-unbalanced trees).",
    //         "%prog [options] [TREE-FILE [TREE-FILE [...]]]");
    // parser.add_option<std::string>(&format, "-f", "--format",
    //         "format for input source ('nexus', 'newick'; default = 'nexus')");
    // parser.add_option<unsigned long>(&num_interpolated_points, "-n", "--profile-size",
    //         "number of interpolated points in profile; if specified, this must be equal or greater to than the largest input data size; if not specified, will default to the largest input data size");
    // parser.add_switch(&calculate_other_distance_metrics,
    //         "-a",
    //         "--calculate-other-distance-metrics",
    //         "calculate other distance metrics (i.e., unlabeled symmetric difference, etc.)");
    // parser.add_option<std::string>(&output_prefix, "-o", "--output-prefix",
    //         "prefix (directory path and filename stem) for output file(s); if not specified, will write to standard output");
    // parser.add_switch(&suppress_header_row, NULL, "--suppress-header-row", "do not write column/field name row in reuslts");
    // // parser.add_switch(&rooted, "--rooted", "force all pairwise_distance_trees to be rooted");
    // // parser.add_switch(&unrooted, "--unrooted", "force all pairwise_distance_trees to be unrooted");
    // parser.add_switch(&quiet, "-q", "--quiet", "suppress all informational/progress messages");
    // parser.parse(argc, argv);

    // // set up logger
    // colugo::Logger logger("pstrudel-trees-refdist");
    // if (quiet) {
    //     logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::WARNING);
    // } else {
    //     logger.add_channel(std::cerr, colugo::Logger::LoggingLevel::INFO);
    // }

    // // get source pairwise_distance_trees
    // std::vector<std::string> args = parser.get_args();
    // std::vector<pstrudel::PairwiseDistanceTree> pairwise_distance_trees;
    // if (args.size() == 0) {
    //     logger.info("(reading trees from standard input)");
    //     pstrudel::treeio::read_from_stream(pairwise_distance_trees, std::cin, format);
    // } else {
    //     unsigned file_idx = 1;
    //     for (auto & arg : args) {
    //         logger.info("Reading tree source file ", file_idx, " of ", args.size(), ": ", arg);
    //         pstrudel::treeio::read_from_filepath(pairwise_distance_trees, arg, format);
    //         ++file_idx;
    //     }
    //     logger.info(pairwise_distance_trees.size(), " trees read.");
    // }

    // // process pairwise_distance_trees
    // logger.info("Calculating tree metrics ...");
    // for (auto & tree : pairwise_distance_trees) {
    //     tree.calc_profile_metrics();
    // }

    // // get number of interpolated points
    // std::map<const std::string, unsigned long> global_num_interpolated_points;
    // if (num_interpolated_points == 0) {
    //     logger.info("Calculating number of interpolated points in profiles ...");
    //     for (auto & tree : pairwise_distance_trees) {
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
    //     for (auto & tree : pairwise_distance_trees) {
    //         tree.set_num_interpolated_profile_points(global_num_interpolated_points);
    //     }
    // }

    // // build profiles
    // logger.info("Building tree profiles ...");
    // for (auto & tree : pairwise_distance_trees) {
    //     tree.build_profile();
    // }

    // ///////////////////////////////////////////////////////////////////////////////////////////////
    // // TODO: handle this more smartly
    // // for now, assume all trees have the same number of leaves
    // // generation of reference trees
    // unsigned long num_leaves = pairwise_distance_trees[0].get_num_leaves();
    // std::vector<pstrudel::PairwiseDistanceNodeValue> leaves;

    // // balanced tree
    // logger.info("Generating balanced tree ..");
    // leaves.reserve(num_leaves);
    // for (unsigned long i = 0; i < num_leaves; ++i) {
    //     leaves.emplace_back();
    // }
    // pstrudel::PairwiseDistanceTree balanced_tree_pwd;
    // build_maximally_balanced_tree(balanced_tree_pwd, leaves.begin(), leaves.end());
    // auto newick_writer1 =  get_newick_writer<pstrudel::PairwiseDistanceTree>(false);
    // std::ostringstream balanced_tree_newick;
    // newick_writer1.write_tree(balanced_tree_pwd, balanced_tree_newick);
    // logger.info("Balanced tree: ", balanced_tree_newick.str());
    // balanced_tree_pwd.calc_profile_metrics();
    // if (num_interpolated_points == 0) {
    //     balanced_tree_pwd.set_num_interpolated_profile_points(global_num_interpolated_points);
    // }
    // balanced_tree_pwd.build_profile();

    // leaves.clear();
    // logger.info("Generating unbalanced tree ..");
    // for (unsigned long i = 0; i < num_leaves; ++i) {
    //     leaves.emplace_back();
    // }
    // pstrudel::PairwiseDistanceTree unbalanced_tree_pwd;
    // build_maximally_unbalanced_tree(unbalanced_tree_pwd, leaves.begin(), leaves.end());
    // std::ostringstream unbalanced_tree_newick;
    // newick_writer1.write_tree(unbalanced_tree_pwd, unbalanced_tree_newick);
    // logger.info("Unbalanced tree: ", unbalanced_tree_newick.str());
    // unbalanced_tree_pwd.calc_profile_metrics();
    // if (num_interpolated_points == 0) {
    //     unbalanced_tree_pwd.set_num_interpolated_profile_points(global_num_interpolated_points);
    // }
    // unbalanced_tree_pwd.build_profile();
    // //
    // ///////////////////////////////////////////////////////////////////////////////////////////////

    // // calculate pairwise profile distances
    // TreeReferenceDistanceType  balanced_tree_profile_distances;
    // TreeReferenceDistanceType  unbalanced_tree_profile_distances;
    // logger.info("Calculating profile distances to reference trees ...");
    // for (unsigned long tree_idx1 = 0; tree_idx1 < pairwise_distance_trees.size(); ++tree_idx1) {
    //     auto & tree1 = pairwise_distance_trees[tree_idx1];
    //     balanced_tree_profile_distances[tree_idx1] = tree1.get_unweighted_subprofile_distance(balanced_tree_pwd);
    //     unbalanced_tree_profile_distances[tree_idx1] = tree1.get_unweighted_subprofile_distance(unbalanced_tree_pwd);
    // }

    // // calculate pairwise symmetric distances
    // TreeReferenceDistanceType  balanced_tree_symdiff_distances;
    // TreeReferenceDistanceType  unbalanced_tree_symdiff_distances;
    // pstrudel::SymmetricDifferenceTree balanced_tree_symd(balanced_tree_pwd);
    // pstrudel::SymmetricDifferenceTree unbalanced_tree_symd(unbalanced_tree_pwd);
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
    //         balanced_tree_symdiff_distances[tree_idx1] = tree1.calc_leaf_set_sizes_unlabeled_symmetric_difference(balanced_tree_symd);
    //         unbalanced_tree_symdiff_distances[tree_idx1] = tree1.calc_leaf_set_sizes_unlabeled_symmetric_difference(unbalanced_tree_symd);
    //     }
    // }

    // std::ostream& out = std::cout;
    // if (!suppress_header_row) {
    //     out << "Tree";
    //     out << "\t" << "profile.dist.to.balanced.tree";
    //     out << "\t" << "profile.dist.to.unbalanced.tree";
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << "sym.diff.dist.to.balanced.tree";
    //         out << "\t" << "sym.diff.dist.to.unbalanced.tree";
    //     }
    //     out << std::endl;
    // }

    // // print distance between reference trees
    // if (!exclude_distance_between_reference_trees) {
    //     out << "-1";
    //     out << "\t" << std::setprecision(20) << balanced_tree_pwd.get_unweighted_subprofile_distance(balanced_tree_pwd);
    //     out << "\t" << std::setprecision(20) << balanced_tree_pwd.get_unweighted_subprofile_distance(unbalanced_tree_pwd);
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << std::setprecision(20) << balanced_tree_symd.calc_leaf_set_sizes_unlabeled_symmetric_difference(balanced_tree_symd);
    //         out << "\t" << std::setprecision(20) << balanced_tree_symd.calc_leaf_set_sizes_unlabeled_symmetric_difference(unbalanced_tree_symd);
    //     }
    //     out << std::endl;
    //     // out << "-2";
    //     // out << "\t" << std::setprecision(20) << unbalanced_tree_pwd.get_unweighted_subprofile_distance(balanced_tree_pwd);
    //     // out << "\t" << std::setprecision(20) << unbalanced_tree_pwd.get_unweighted_subprofile_distance(unbalanced_tree_pwd);
    //     // out << std::endl;
    // }

    // unsigned long num_trees = pairwise_distance_trees.size();
    // for (unsigned long tree_idx1 = 0; tree_idx1 < num_trees; ++tree_idx1) {
    //     out << tree_idx1 + 1;
    //     out << "\t" << std::setprecision(20) << balanced_tree_profile_distances[tree_idx1];
    //     out << "\t" << std::setprecision(20) << unbalanced_tree_profile_distances[tree_idx1];
    //     if (calculate_other_distance_metrics) {
    //         out << "\t" << balanced_tree_symdiff_distances[tree_idx1];
    //         out << "\t" << unbalanced_tree_symdiff_distances[tree_idx1];
    //     }
    //     out << std::endl;
    // }

}


