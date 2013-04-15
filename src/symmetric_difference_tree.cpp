#include "symmetric_difference_tree.hpp"

///////////////////////////////////////////////////////////////////////////////
// SymmetricDifferenceTree

namespace pstrudel {

SymmetricDifferenceTree::SymmetricDifferenceTree() {
}

SymmetricDifferenceTree::SymmetricDifferenceTree(const SymmetricDifferenceTree & other)
    : BasicTree<SymmetricDifferenceNodeValue>(other)
      , subtree_leaf_set_sizes_(other.subtree_leaf_set_sizes_)
      , subtree_clade_sizes_(other.subtree_clade_sizes_) {
}

void SymmetricDifferenceTree::calc_subtree_sizes() {
    this->subtree_leaf_set_sizes_.clear();
    for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        if (ndi.is_leaf()) {
            ndi->set_num_leaves(1);
            ndi->set_num_descendent_nodes(0);
        } else {
            unsigned long nleaves = 0;
            unsigned long ndesc = 0;
            for (auto chi = this->children_begin(ndi); chi != this->children_end(ndi); ++chi) {
                nleaves += chi->get_num_leaves();
                ndesc += chi->get_num_descendent_nodes();
                ndesc += 1;
                // std::cerr << nleaves << std::endl;
            }
            ndi->set_num_leaves(nleaves);
            ndi->set_num_descendent_nodes(ndesc);
            this->subtree_leaf_set_sizes_.insert(nleaves);
            this->subtree_clade_sizes_.insert(ndesc);
        }
    }
}

unsigned long SymmetricDifferenceTree::calc_leaf_set_sizes_unlabeled_symmetric_difference(SymmetricDifferenceTree & other) {
    if (this->subtree_leaf_set_sizes_.empty()) {
        this->calc_subtree_sizes();
    }
    if (other.subtree_leaf_set_sizes_.empty()) {
        other.calc_subtree_sizes();
    }
    return SymmetricDifferenceTree::calc_set_symmetric_difference(this->subtree_leaf_set_sizes_, other.subtree_leaf_set_sizes_);
}

unsigned long SymmetricDifferenceTree::calc_clade_sizes_unlabeled_symmetric_difference(SymmetricDifferenceTree & other) {
    if (this->subtree_clade_sizes_.empty()) {
        this->calc_subtree_sizes();
    }
    if (other.subtree_clade_sizes_.empty()) {
        other.calc_subtree_sizes();
    }
    return SymmetricDifferenceTree::calc_set_symmetric_difference(this->subtree_clade_sizes_, other.subtree_clade_sizes_);
}

unsigned long SymmetricDifferenceTree::get_unlabeled_symmetric_difference(SymmetricDifferenceTree & other) {
    return this->calc_leaf_set_sizes_unlabeled_symmetric_difference(other);
}

unsigned long SymmetricDifferenceTree::get_unweighted_labeled_symmetric_difference(SymmetricDifferenceTree & other) {
    throw std::logic_error("Not Implemented: SymmetricDifferenceTree::get_unweighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long SymmetricDifferenceTree::get_weighted_labeled_symmetric_difference(SymmetricDifferenceTree & other) {
    throw std::logic_error("Not Implemented: SymmetricDifferenceTree::get_weighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long SymmetricDifferenceTree::calc_set_symmetric_difference(
        const SizesSetType & set1,
        const SizesSetType & set2,
        SizesSetType * common,
        SizesSetType * unmatched1,
        SizesSetType * unmatched2) {
    auto s1 = set1;
    auto s2 = set2;
    auto s1_iter = s1.begin();
    while (!s1.empty() && !s2.empty() && s1_iter != s1.end()) {
        auto s2_iter = s2.find(*s1_iter);
        if (s2_iter != s2.end()) {
            if (common != nullptr) {
                common->insert(*s1_iter);
            }
            s1.erase(s1_iter++);
            s2.erase(s2_iter);
        } else {
            ++s1_iter;
        }
    }
    if (unmatched1 !=  nullptr) {
        unmatched1->insert(s1.begin(), s1.end());
    }
    if (unmatched2 !=  nullptr) {
        unmatched2->insert(s2.begin(), s2.end());
    }
    return s1.size() + s2.size();
}

}
