#ifndef PSTRUDEL_DATAIO_HPP
#define PSTRUDEL_DATAIO_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <functional>
#include <map>
#include <ncl/nxsmultiformat.h>
#include <colugo/console.hpp>
#include <platypus/parse/nclreader.hpp>
#include <platypus/model/standardinterface.hpp>

namespace pstrudel {

namespace sequenceio {

////////////////////////////////////////////////////////////////////////////////
// Sequence Reading Utility Functions

template <class SequencesType>
int read_from_stream(
        SequencesType& sequences,
        std::istream& src,
        const std::string& format) {
    MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
    reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
    reader.SetCoerceUnderscoresToSpaces(true);
    const char * format_cstr = nullptr;
    if (format == "phylip" || format == "dnaphylip") {
        format_cstr = "dnarelaxedphylip";
    } else if (format == "rnaphylip") {
        format_cstr = "rnarelaxedphylip";
    } else if (format == "fasta") {
        format_cstr = "dnafasta";
    } else {
        format_cstr = format.c_str();
    }
    reader.ReadStream(src, format_cstr);
    unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
    NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
    if (!taxa_block) {
        colugo::console::abort("No taxon definitions were parsed (invalid file format?)");
    }
    NxsCharactersBlock * chars_block = reader.GetCharactersBlock(taxa_block, 0);
    if (!chars_block) {
        colugo::console::abort("No character states were parsed (invalid file format?)");
    }
    unsigned int ntax = taxa_block->GetNTaxTotal();
    for (unsigned int taxon_idx = 0; taxon_idx < ntax; ++taxon_idx) {
        // const char * label = NxsString::GetEscaped(taxa_block->GetTaxonLabel(taxon_idx)).c_str();
        const char * label = taxa_block->GetTaxonLabel(taxon_idx).c_str();
        NxsDiscreteStateRow row = chars_block->GetDiscreteMatrixRow(taxon_idx);
        // seqs->reserve(row.size());
        std::ostringstream o;
        for (unsigned int col = 0; col < row.size(); ++col) {
            chars_block->ShowStateLabels(o, taxon_idx, col, 0);
        }
        auto * seq = sequences.new_sequence(label);
        seq->append_states_by_symbols(o.str());
        // for (auto & ncl_state : row) {
        //     seqs->append_state_by_symbol(symbols[ncl_state]);
        // }
    }
    return ntax;
}

template <class SequencesType>
int read_from_filepath(
        SequencesType& sequences,
        const std::string& filepath,
        const std::string& format) {
    std::ifstream f(filepath);
    if (!f.good()) {
        colugo::console::abort("Error opening file for input");
    }
    return read_from_stream(sequences, f, format);
}

template <class SequencesType>
int read_from_string(
        SequencesType& sequences,
        const std::string& str,
        const std::string& format) {
    std::istringstream s(str);
    return read_from_stream(sequences, s, format);
}

////////////////////////////////////////////////////////////////////////////////
// Sequence Writing Utility Functions

template <class SequencesType>
int write_fasta(
        SequencesType& sequences,
        std::ostream& out) {
    // PSTRUDEL_ASSERT(this->sequences_.size() == this->labels_.size());
    // for (unsigned long idx = 0; idx < this->sequences_.size(); ++idx) {
    //     out << ">" << this->labels_[idx] << std::endl;
    //     unsigned int col_count = 0;
    //     for (auto & ch : (*this->sequences_[idx]) ) {
    //         if (col_count == 70) {
    //             out << std::endl;
    //             col_count = 0;
    //         }
    //         out << NucleotideSequences::get_symbol_from_state(ch);
    //         ++col_count;
    //     }
    //     out << "\n" << std::endl;
    // }
    return 1;
}

} // namespace sequenceio

namespace treeio {

// convenience functions

template <class TreeT>
void postprocess_tree(TreeT & tree, unsigned long idx, unsigned long ntips, unsigned long nints, double tree_length) {
    // tree.set_index(idx);
    tree.set_num_tips(ntips);
    // tree.set_nints(nints);
    tree.set_total_tree_length(tree_length);
}

template <class TreeT>
int read_from_stream(std::vector<TreeT>& trees,
        std::istream& src,
        const std::string& format) {
    auto reader = platypus::NclTreeReader<TreeT>();
    platypus::bind_standard_interface(reader);
    reader.set_tree_postprocess_fn(postprocess_tree<TreeT>);
    std::function<TreeT& ()> get_new_tree_reference = [&trees] () -> TreeT& { trees.emplace_back(); return trees.back(); };
    return reader.read(src, get_new_tree_reference);
}

template <class TreeT>
int read_from_filepath(std::vector<TreeT>& trees,
        const std::string& filepath,
        const std::string& format) {
    std::ifstream f(filepath);
    if (!f.good()) {
        throw platypus::ReaderException(__FILE__, __LINE__, "platypus::BaseTreeReader::read_from_filepath(): Error opening file for input");
    }
    return read_from_stream(trees, f, format);
}

template <class TreeT>
int read_from_string(std::vector<TreeT>& trees,
        const std::string& str,
        const std::string& format) {
    std::istringstream s(str);
    return read_from_stream(trees, s, format);
}


////////////////////////////////////////////////////////////////////////////////
// Tree Writing Utility Functions

template <class TreeT, class iter>
void write_newick_node(const TreeT& tree,
        const iter& tree_iter,
        std::ostream& out,
        bool include_edge_lens=true,
        unsigned long edge_len_precision=10) {
    COLUGO_NDEBUG_ASSERT(tree_iter);
    if (!tree_iter.is_leaf()) {
        out << "(";
        int ch_count = 0;
        for (auto chi = tree.children_begin(tree_iter);
                chi != tree.children_end(tree_iter);
                ++chi, ++ch_count) {
            if (ch_count > 0) {
                out << ", ";
            }
            write_newick_node(tree, chi, out, include_edge_lens, edge_len_precision);
        }
        out << ")";
    } else {
    }
    // label = NxsString::GetEscaped(tb->GetTaxonLabel(ncl_taxon_idx)).c_str();
    // out << tree_iter->get_label();
    out << NxsString::GetEscaped(tree_iter->get_label()).c_str();
    if (include_edge_lens) {
        out << ":" << std::setprecision(edge_len_precision) << tree_iter->get_edge_length();
    }
}

template <class TreeT>
void write_newick(const TreeT& tree,
        std::ostream& out,
        bool include_edge_lens=true,
        unsigned long edge_len_precision=10) {
    if (tree.is_rooted()) {
        out << "[&R] ";
    } else {
        out << "[&U] ";
    }
    write_newick_node(tree, tree.begin(), out, include_edge_lens, edge_len_precision);
    out << ";" << std::endl;
}

template <class TreeT>
void write_nexus_tree(const TreeT& tree,
        std::ostream& out,
        bool include_edge_lens=true,
        unsigned long edge_len_precision=10) {
    out << "    TREE ";
    if (tree.get_label().empty()) {
        out << "unlabeled";
    } else {
        out << tree.get_label();
    }
    out << " = ";
    if (tree.is_rooted()) {
        out << "[&R] ";
    } else {
        out << "[&U] ";
    }
    write_newick_node(tree, tree.begin(), out, include_edge_lens, edge_len_precision);
    out << ";\n";
}

template <class LabelContainer>
void write_nexus_taxa_block(const LabelContainer& tax_labels, std::ostream& out) {
    out << "BEGIN TAXA;\n";
    out << "    DIMENSIONS NTAX=" << tax_labels.size() << ";\n";
    out << "    TAXLABELS\n";
    for (auto & tax_label : tax_labels) {
        out << "        " << NxsString::GetEscaped(tax_label.c_str()) << "\n";
    }
    out << "    ;\n";
    out << "END;\n\n";
}

template <class TreeT>
void write_nexus(std::vector<TreeT>& trees,
        std::ostream& out,
        bool include_edge_lens=true,
        unsigned long edge_len_precision=10) {
    out << "#NEXUS\n\n";
    if (trees.size() == 0) {
        return;
    }
    std::set<std::string> tax_labels;
    auto & first_tree = trees.at(0);
    for (auto ndi = first_tree.leaf_begin(); ndi != first_tree.leaf_end(); ++ndi) {
        tax_labels.insert(ndi->get_label());
    }
    write_nexus_taxa_block(tax_labels, out);
    out << "BEGIN TREES;\n";
    for (auto & tree : trees) {
        write_nexus_tree(tree, out, include_edge_lens, edge_len_precision);
    }
    out << "END;\n\n";
}

} // namespace treeio

} // namespace pstrudel

#endif

// // template <class TreeT, class NodeLabelType, class EdgeLenType>
// // int read_from_stream(
// //         std::istream& src,
// //         const std::string& format,
// //         std::function<TreeT& ()> tree_factory,
// //         void (TreeT::value_type::*node_data_label_setter)(const NodeLabelType&)=nullptr,
// //         void (TreeT::value_type::*node_data_edge_setter)(EdgeLenType)=nullptr,
// //         void (TreeT::*tree_is_rooted_setter)(bool)=nullptr,
// //         typename TreeT::node_type * (TreeT::*leaf_node_factory)()=nullptr,
// //         typename TreeT::node_type * (TreeT::*internal_node_factory)()=nullptr,
// //         unsigned long max_tips=0) {
// //     MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
// //     reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
// //     reader.SetCoerceUnderscoresToSpaces(false);
// //     const char * format_cstr = nullptr;
// //     if (format == "newick") {
// //         format_cstr = "relaxedphyliptree";
// //     } else {
// //         format_cstr = format.c_str();
// //     }
// //     reader.ReadStream(src, format_cstr);
// //     unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
// //     NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
// //     if (!taxa_block) {
// //         colugo::console::abort("No taxon definitions were parsed (invalid file format?)");
// //     }
// //     NxsTreesBlock * trees_block = reader.GetTreesBlock(taxa_block, 0);
// //     if (!trees_block) {
// //         return 0;
// //     }
// //     if (max_tips == 0) {
// //         max_tips = taxa_block->GetNTaxTotal();
// //     }
// //     unsigned int num_trees = trees_block->GetNumTrees();
// //     int tree_count = 0;
// //     for (unsigned int tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
// //         auto & tree = tree_factory();
// //         const NxsFullTreeDescription & ftd = trees_block->GetFullTreeDescription(tree_idx);
// //         build_tree<TreeT>(tree,
// //                 taxa_block,
// //                 ftd,
// //                 node_data_label_setter,
// //                 node_data_edge_setter,
// //                 tree_is_rooted_setter,
// //                 leaf_node_factory,
// //                 internal_node_factory);
// //         ++tree_count;
// //     }
// //     return tree_count;
// // }

// pstrudel::treeio::read_from_stream<pstrudel::PairwiseDistanceTree>(
//         src,
//         format,
//         [&trees] () -> pstrudel::PairwiseDistanceTree& { trees.emplace_back(); return trees.back(); },
//         &pstrudel::PairwiseDistanceTree::value_type::set_label,
//         &pstrudel::PairwiseDistanceTree::value_type::set_edge_length);

// pstrudel::treeio::read_from_stream<pstrudel::PairwiseDistanceTree>(
//         src,
//         format,
//         [&trees] () -> pstrudel::PairwiseDistanceTree& { trees.emplace_back(); return trees.back(); },
//         [] (pstrudel::PairwiseDistanceTree& tree, bool is_rooted) { tree.set_rooted(is_rooted); },
//         [] (pstrudel::PairwiseDistanceTree::value_type& node, const std::string& label) { node.set_label(label); },
//         [] (pstrudel::PairwiseDistanceTree::value_type& node, double edge_length) { node.set_edge_length(edge_length); },
//         [] (pstrudel::PairwiseDistanceTree& tree) -> pstrudel::PairwiseDistanceTree::node_type * {return tree.create_node(); },
//         [] (pstrudel::PairwiseDistanceTree& tree) -> pstrudel::PairwiseDistanceTree::node_type * {return tree.create_node(); }
//         );

// pstrudel::treeio::read_from_stream<pstrudel::PairwiseDistanceTree>(
//         src,
//         format,
//         [&trees] () -> pstrudel::PairwiseDistanceTree& { trees.emplace_back(); return trees.back(); },
//         [] (pstrudel::PairwiseDistanceTree& tree, bool is_rooted) { tree.set_rooted(is_rooted); },
//         [] (pstrudel::PairwiseDistanceTree::value_type& node, const std::string& label) { node.set_label(label); },
//         [] (pstrudel::PairwiseDistanceTree::value_type& node, double edge_length) { node.set_edge_length(edge_length); }
//         );

/**
 * Dummy function to serve as default template parameter in operational
 * functions below.
 */
// class TreenaryDummySetter_ {
//     public:
//         template <typename... Types>
//         void operator()(const Types&... args) {
//         }
// };

/**
 * Default node allocator.
 */
// class TreenaryDefaultNodeFactory_ {
//     public:
//         template <typename TreeT>
//         auto operator()(TreeT& tree) -> decltype(tree.create_node()) {
//             return tree.create_node();
//         }
// };

// template <class TreeT,
//         class IsRootedFuncPtrType=TreenaryDummySetter_,
//         class NodeLabelFuncPtrType=TreenaryDummySetter_,
//         class NodeEdgeFuncPtrType=TreenaryDummySetter_,
//         class LeafNodeFactoryType=TreenaryDefaultNodeFactory_,
//         class InternalNodeFactoryType=TreenaryDefaultNodeFactory_>
// class TreeReader {

//     public:
//         TreeReader(
//                 std::function<TreeT& ()> tree_factory,
//                 IsRootedFuncPtrType tree_is_rooted_setter=IsRootedFuncPtrType(),
//                 NodeLabelFuncPtrType node_data_label_setter=NodeLabelFuncPtrType(),
//                 NodeEdgeFuncPtrType node_data_edge_setter=NodeEdgeFuncPtrType(),
//                 LeafNodeFactoryType leaf_node_factory=LeafNodeFactoryType(),
//                 InternalNodeFactoryType internal_node_factory=InternalNodeFactoryType())
//                 : tree_factory_(tree_factory)
//                   , tree_is_rooted_setter_(tree_is_rooted_setter)
//                   , node_data_label_setter_(node_data_label_setter)
//                   , node_data_edge_setter_(node_data_label_setter)
//                   , leaf_node_factory_(leaf_node_factory_)
//                   , internal_node_factory_(internal_node_factory_) {
//         }

//         void set_node_value_label_setter(NodeLabelFuncPtrType setter) {
//             this->node_data_label_setter_ = setter;
//         }

//     private:
//         std::function<TreeT& ()>     tree_factory_;
//         IsRootedFuncPtrType              tree_is_rooted_setter_;
//         NodeLabelFuncPtrType             node_data_label_setter_;
//         NodeEdgeFuncPtrType              node_data_edge_setter_;
//         LeafNodeFactoryType             leaf_node_factory_;
//         InternalNodeFactoryType         internal_node_factory_;

// }; // TreeReader


// template <class TreeT,
//         class IsRootedFuncPtrType=TreenaryDummySetter_,
//         class NodeLabelFuncPtrType=TreenaryDummySetter_,
//         class NodeEdgeFuncPtrType=TreenaryDummySetter_,
//         class LeafNodeFactoryType=TreenaryDefaultNodeFactory_,
//         class InternalNodeFactoryType=TreenaryDefaultNodeFactory_>
// int build_tree(
//         TreeT& ttree,
//         const NxsTaxaBlock * tb,
//         const NxsFullTreeDescription & ftd,
//         IsRootedFuncPtrType is_rooted_setter,
//         NodeLabelFuncPtrType node_label_setter,
//         NodeEdgeFuncPtrType node_edge_setter,
//         LeafNodeFactoryType leaf_node_factory,
//         InternalNodeFactoryType internal_node_factory) {
//     is_rooted_setter(ttree, ftd.IsRooted());
//     NxsSimpleTree ncl_tree(ftd, -1, -1.0);
//     auto * root = ttree.head_node();
//     decltype(root) node_parent = nullptr;
//     decltype(root) new_node = nullptr;
//     std::vector<const NxsSimpleNode *> ncl_nodes = ncl_tree.GetPreorderTraversal();
//     std::map<const NxsSimpleNode *, decltype(root)> ncl_to_native;
//     int size = 0;
//     for (auto & ncl_node : ncl_nodes) {
//         const NxsSimpleEdge & ncl_edge = ncl_node->GetEdgeToParentRef();
//         const NxsSimpleNode * ncl_par = ncl_edge.GetParent();
//         std::vector<NxsSimpleNode *> ncl_child_nodes = ncl_node->GetChildren();
//         unsigned int nchildren = ncl_child_nodes.size();
//         double edge_len = ncl_edge.GetDblEdgeLen();
//         if (edge_len < 0) {
//             edge_len = 0.0;
//         }
//         std::string label;

//         // if (nchildren > 2) {
//         //     colugo::console::abort("Tree has node with more than 2 children");
//         if (nchildren == 1) {
//             colugo::console::abort("Tree has node with only 1 child");
//         } else if (nchildren == 0) {
//             unsigned int ncl_taxon_idx = ncl_node->GetTaxonIndex();
//             // label = NxsString::GetEscaped(tb->GetTaxonLabel(ncl_taxon_idx)).c_str();
//             label = tb->GetTaxonLabel(ncl_taxon_idx).c_str();
//             // new_node = ttree.allocate_leaf_node();
//             // ((obj).*(f1))();
//             new_node = leaf_node_factory(ttree);
//         } else {
//             label = NxsString::GetEscaped(ncl_node->GetName()).c_str();
//             if (!ncl_par) {
//                 new_node = ttree.head_node();
//             } else {
//                 // new_node = ttree.allocate_internal_node();
//                 new_node = internal_node_factory(ttree);
//             }
//         }
//         node_label_setter(new_node->data(), label);
//         node_edge_setter(new_node->data(), edge_len);
//         ncl_to_native[ncl_node] = new_node;
//         if (ncl_par) {
//             if (ncl_to_native.find(ncl_par) == ncl_to_native.end()) {
//                 colugo::console::abort("Parent node not visited in preorder traversal");
//             }
//             node_parent = ncl_to_native[ncl_par];
//             if (!node_parent) {
//                 colugo::console::abort("Null parent node");
//             }
//             node_parent->add_child(new_node);
//         }
//         ++size;
//     }
//     return size;
// }

// template <class TreeT,
//         class IsRootedFuncPtrType=TreenaryDummySetter_,
//         class NodeLabelFuncPtrType=TreenaryDummySetter_,
//         class NodeEdgeFuncPtrType=TreenaryDummySetter_,
//         class LeafNodeFactoryType=TreenaryDefaultNodeFactory_,
//         class InternalNodeFactoryType=TreenaryDefaultNodeFactory_>
// int read_from_stream(
//         std::istream& src,
//         const std::string& format,
//         std::function<TreeT& ()> tree_factory,
//         IsRootedFuncPtrType tree_is_rooted_setter=IsRootedFuncPtrType(),
//         NodeLabelFuncPtrType node_data_label_setter=NodeLabelFuncPtrType(),
//         NodeEdgeFuncPtrType node_data_edge_setter=NodeEdgeFuncPtrType(),
//         LeafNodeFactoryType leaf_node_factory=LeafNodeFactoryType(),
//         InternalNodeFactoryType internal_node_factory=InternalNodeFactoryType()) {
//     MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
//     reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
//     reader.SetCoerceUnderscoresToSpaces(false);
//     const char * format_cstr = nullptr;
//     if (format == "newick") {
//         format_cstr = "relaxedphyliptree";
//     } else {
//         format_cstr = format.c_str();
//     }
//     reader.ReadStream(src, format_cstr);
//     unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
//     NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
//     if (!taxa_block) {
//         colugo::console::abort("No taxon definitions were parsed (invalid file format?)");
//     }
//     NxsTreesBlock * trees_block = reader.GetTreesBlock(taxa_block, 0);
//     if (!trees_block) {
//         return 0;
//     }
//     unsigned int num_trees = trees_block->GetNumTrees();
//     int tree_count = 0;
//     for (unsigned int tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
//         auto & tree = tree_factory();
//         const NxsFullTreeDescription & ftd = trees_block->GetFullTreeDescription(tree_idx);
//         build_tree(tree,
//                 taxa_block,
//                 ftd,
//                 tree_is_rooted_setter,
//                 node_data_label_setter,
//                 node_data_edge_setter,
//                 leaf_node_factory,
//                 internal_node_factory);
//         ++tree_count;
//     }
//     return tree_count;
// }



/**
 * Constructs a ternary-derived tree, ``TreeT``, from NCL data objects.
 *
 * @param  ttree                    A ternary-derived Tree object, assumed to
 *                                  be initialized (i.e., with head node and
 *                                  stop node allocated).
 *
 * @param  node_data_label_setter   Pointer to member function of a the node
 *                                  data class (``TreeT::value_type``) that
 *                                  accepts a std::string reference argument
 *                                  representing the label of the node (or
 *                                  taxon, if a leaf node) and sets the
 *                                  ``TreeT::value_type`` object state
 *                                  accordingly.
 *
 * @param  node_data_edge_setter    Pointer to member function of a the node
 *                                  data class (``TreeT::value_type``) that
 *                                  accepts a double argument representing the
 *                                  length of the edge
 *                                  subtending the node and sets the
 *                                  ``TreeT::value_type`` object state
 *                                  accordingly.
 *
 * @param  tree_is_rooted           Pointer to member function of a the
 *                                  TreeT class that accepts a boolean
 *                                  argument representing the rooted state of
 *                                  the tree (``true`` if the tree is rooted,
 *                                  ``false`` otherwise) and sets the state of
 *                                  the TreeT object accordingly.
 *
 * @param  leaf_node_factory        Pointer to member function of a the
 *                                  TreeT class that allocates, constructs
 *                                  and returns a pointer to a new leaf node.
 *                                  If not specified, then the default
 *                                  platypus allocation and construction
 *                                  function, ``Tree::create_node()`` will be used.
 *
 * @param  leaf_node_factory        Pointer to member function of a the
 *                                  TreeT class that allocates, constructs
 *                                  and returns a pointer to a new internal
 *                                  node.  If not specified, then the default
 *                                  platypus allocator will be used.
 *
 * @returns                         The number of nodes on the tree.
 */
// template <class TreeT, class NodeLabelType, class EdgeLenType>
// int build_tree(
//         TreeT& ttree,
//         const NxsTaxaBlock * tb,
//         const NxsFullTreeDescription & ftd,
//         void (TreeT::value_type::*node_data_label_setter)(const NodeLabelType&)=nullptr,
//         void (TreeT::value_type::*node_data_edge_setter)(EdgeLenType)=nullptr,
//         void (TreeT::*tree_is_rooted_setter)(bool)=nullptr,
//         typename TreeT::node_type * (TreeT::*leaf_node_factory)()=nullptr,
//         typename TreeT::node_type * (TreeT::*internal_node_factory)()=nullptr) {
//     if (tree_is_rooted_setter != nullptr) {
//         ((ttree).*(tree_is_rooted_setter))(ftd.IsRooted());
//     }
//     if (leaf_node_factory == nullptr) {
//         leaf_node_factory = &TreeT::create_node;
//     }
//     if (internal_node_factory == nullptr) {
//         internal_node_factory = &TreeT::create_node;
//     }
//     NxsSimpleTree ncl_tree(ftd, -1, -1.0);
//     auto * root = ttree.head_node();
//     decltype(root) node_parent = nullptr;
//     decltype(root) new_node = nullptr;
//     std::vector<const NxsSimpleNode *> ncl_nodes = ncl_tree.GetPreorderTraversal();
//     std::map<const NxsSimpleNode *, decltype(root)> ncl_to_native;
//     int size = 0;
//     for (auto & ncl_node : ncl_nodes) {
//         const NxsSimpleEdge & ncl_edge = ncl_node->GetEdgeToParentRef();
//         const NxsSimpleNode * ncl_par = ncl_edge.GetParent();
//         std::vector<NxsSimpleNode *> ncl_child_nodes = ncl_node->GetChildren();
//         unsigned int nchildren = ncl_child_nodes.size();
//         double edge_len = ncl_edge.GetDblEdgeLen();
//         if (edge_len < 0) {
//             edge_len = 0.0;
//         }
//         std::string label;

//         // if (nchildren > 2) {
//         //     colugo::console::abort("Tree has node with more than 2 children");
//         if (nchildren == 1) {
//             colugo::console::abort("Tree has node with only 1 child");
//         } else if (nchildren == 0) {
//             unsigned int ncl_taxon_idx = ncl_node->GetTaxonIndex();
//             // label = NxsString::GetEscaped(tb->GetTaxonLabel(ncl_taxon_idx)).c_str();
//             label = tb->GetTaxonLabel(ncl_taxon_idx).c_str();
//             // new_node = ttree.allocate_leaf_node();
//             // ((obj).*(f1))();
//             new_node = ((ttree).*(leaf_node_factory))();
//         } else {
//             label = NxsString::GetEscaped(ncl_node->GetName()).c_str();
//             if (!ncl_par) {
//                 new_node = ttree.head_node();
//             } else {
//                 // new_node = ttree.allocate_internal_node();
//                 new_node = ((ttree).*(internal_node_factory))();
//             }
//         }
//         if (node_data_label_setter != nullptr) {
//             ((new_node->data()).*(node_data_label_setter))(label);
//         }
//         if (node_data_edge_setter != nullptr) {
//             ((new_node->data()).*(node_data_edge_setter))(edge_len);
//         }
//         ncl_to_native[ncl_node] = new_node;
//         if (ncl_par) {
//             if (ncl_to_native.find(ncl_par) == ncl_to_native.end()) {
//                 colugo::console::abort("Parent node not visited in preorder traversal");
//             }
//             node_parent = ncl_to_native[ncl_par];
//             if (!node_parent) {
//                 colugo::console::abort("Null parent node");
//             }
//             node_parent->add_child(new_node);
//         }
//         ++size;
//     }
//     return size;
// }

// template <class TreeT, class NodeLabelType, class EdgeLenType>
// int read_from_stream(
//         std::vector<TreeT>& trees,
//         std::istream& src,
//         const std::string& format,
//         void (TreeT::value_type::*node_data_label_setter)(const NodeLabelType&)=nullptr,
//         void (TreeT::value_type::*node_data_edge_setter)(EdgeLenType)=nullptr,
//         void (TreeT::*tree_is_rooted_setter)(bool)=nullptr,
//         typename TreeT::node_type * (TreeT::*leaf_node_factory)()=nullptr,
//         typename TreeT::node_type * (TreeT::*internal_node_factory)()=nullptr,
//         unsigned long max_tips=0) {
//     MultiFormatReader reader(-1, NxsReader::IGNORE_WARNINGS);
//     reader.SetWarningOutputLevel(NxsReader::AMBIGUOUS_CONTENT_WARNING);
//     reader.SetCoerceUnderscoresToSpaces(false);
//     const char * format_cstr = nullptr;
//     if (format == "newick") {
//         format_cstr = "relaxedphyliptree";
//     } else {
//         format_cstr = format.c_str();
//     }
//     reader.ReadStream(src, format_cstr);
//     unsigned num_taxa_blocks = reader.GetNumTaxaBlocks();
//     NxsTaxaBlock *  taxa_block = reader.GetTaxaBlock(num_taxa_blocks-1);
//     if (!taxa_block) {
//         colugo::console::abort("No taxon definitions were parsed (invalid file format?)");
//     }
//     NxsTreesBlock * trees_block = reader.GetTreesBlock(taxa_block, 0);
//     if (!trees_block) {
//         return 0;
//     }
//     if (max_tips == 0) {
//         max_tips = taxa_block->GetNTaxTotal();
//     }
//     unsigned int num_trees = trees_block->GetNumTrees();
//     int tree_count = 0;
//     for (unsigned int tree_idx = 0; tree_idx < num_trees; ++tree_idx) {
//         trees.emplace_back(max_tips);
//         auto & tree = trees.back();
//         const NxsFullTreeDescription & ftd = trees_block->GetFullTreeDescription(tree_idx);
//         build_tree<TreeT>(tree,
//                 taxa_block,
//                 ftd,
//                 node_data_label_setter,
//                 node_data_edge_setter,
//                 tree_is_rooted_setter,
//                 leaf_node_factory,
//                 internal_node_factory);
//         ++tree_count;
//     }
//     return tree_count;
// }

// template <class TreeT, class NodeLabelType, class EdgeLenType>
// int read_from_filepath(
//         std::vector<TreeT>& trees,
//         const std::string& filepath,
//         const std::string& format,
//         void (TreeT::value_type::*node_data_label_setter)(const NodeLabelType&)=nullptr,
//         void (TreeT::value_type::*node_data_edge_setter)(EdgeLenType)=nullptr,
//         void (TreeT::*tree_is_rooted_setter)(bool)=nullptr,
//         typename TreeT::node_type * (TreeT::*leaf_node_factory)()=nullptr,
//         typename TreeT::node_type * (TreeT::*internal_node_factory)()=nullptr,
//         unsigned long max_tips=0) {
//     std::ifstream f(filepath);
//     if (!f.good()) {
//         colugo::console::abort("Error opening file for input");
//     }
//     return read_from_stream<TreeT>(trees,
//             f,
//             format,
//             node_data_label_setter,
//             node_data_edge_setter,
//             tree_is_rooted_setter,
//             leaf_node_factory,
//             internal_node_factory,
//             max_tips);
// }

// template <class TreeT, class NodeLabelType, class EdgeLenType>
// int read_from_string(
//         std::vector<TreeT>& trees,
//         const std::string& str,
//         const std::string& format,
//         void (TreeT::value_type::*node_data_label_setter)(const NodeLabelType&)=nullptr,
//         void (TreeT::value_type::*node_data_edge_setter)(EdgeLenType)=nullptr,
//         void (TreeT::*tree_is_rooted_setter)(bool)=nullptr,
//         typename TreeT::node_type * (TreeT::*leaf_node_factory)()=nullptr,
//         typename TreeT::node_type * (TreeT::*internal_node_factory)()=nullptr,
//         unsigned long max_tips=0) {
//     std::istringstream s(str);
//     return read_from_stream<TreeT>(trees,
//             s,
//             format,
//             node_data_label_setter,
//             node_data_edge_setter,
//             tree_is_rooted_setter,
//             leaf_node_factory,
//             internal_node_factory,
//             max_tips);
// }

