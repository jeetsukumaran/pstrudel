#include "distancetree.hpp"

namespace pstrudel {

//////////////////////////////////////////////////////////////////////////////
// DistanceTreeProfile


//////////////////////////////////////////////////////////////////////////////
// DistanceTreeSymmetricDifference

//////////////////////////////////////////////////////////////////////////////
// DistanceTree

DistanceTree::DistanceTree(bool is_rooted)
        : platypus::StandardTree<DistanceNodeValue>(is_rooted) {
}

DistanceTree::DistanceTree(DistanceTree && other)
    : platypus::StandardTree<DistanceNodeValue>(other)
    , unweighted_pairwise_tip_distance_profile_(std::move(other.unweighted_pairwise_tip_distance_profile_))
    , weighted_pairwise_tip_distance_profile_(std::move(other.weighted_pairwise_tip_distance_profile_))
    , number_of_tips_(other.number_of_tips_)
    , total_tree_length_(other.total_tree_length_)
    , subtree_leaf_set_sizes_(std::move(other.subtree_leaf_set_sizes_)) {
}

DistanceTree::DistanceTree(const DistanceTree & other)
    : platypus::StandardTree<DistanceNodeValue>(other)
    , unweighted_pairwise_tip_distance_profile_(other.unweighted_pairwise_tip_distance_profile_)
    , weighted_pairwise_tip_distance_profile_(other.weighted_pairwise_tip_distance_profile_)
    , number_of_tips_(other.number_of_tips_)
    , total_tree_length_(other.total_tree_length_)
    , subtree_leaf_set_sizes_(other.subtree_leaf_set_sizes_) {
}

DistanceTree & DistanceTree::operator=(const DistanceTree & other) {
    platypus::StandardTree<DistanceNodeValue>::operator=(other);
    unweighted_pairwise_tip_distance_profile_ = other.unweighted_pairwise_tip_distance_profile_;
    weighted_pairwise_tip_distance_profile_ = other.weighted_pairwise_tip_distance_profile_;
    this->number_of_tips_ = other.number_of_tips_;
    this->total_tree_length_ = other.total_tree_length_;
    this->subtree_leaf_set_sizes_ = other.subtree_leaf_set_sizes_;
    return *this;
}

void DistanceTree::calc_subtree_sizes() {
    this->subtree_leaf_set_sizes_.clear();
    for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        if (ndi.is_leaf()) {
            ndi->set_num_leaves(1);
        } else {
            unsigned long nleaves = 0;
            unsigned long ndesc = 0;
            for (auto chi = this->children_begin(ndi); chi != this->children_end(ndi); ++chi) {
                nleaves += chi->get_num_leaves();
                ndesc += 1;
            }
            ndi->set_num_leaves(nleaves);
            this->subtree_leaf_set_sizes_.insert(nleaves);
        }
    }
}

unsigned long DistanceTree::calc_leaf_set_sizes_unlabeled_symmetric_difference(DistanceTree & other) {
    if (this->subtree_leaf_set_sizes_.empty()) {
        this->calc_subtree_sizes();
    }
    if (other.subtree_leaf_set_sizes_.empty()) {
        other.calc_subtree_sizes();
    }
    return DistanceTree::calc_set_symmetric_difference(this->subtree_leaf_set_sizes_, other.subtree_leaf_set_sizes_);
}

unsigned long DistanceTree::get_unlabeled_symmetric_difference(DistanceTree & other) {
    return this->calc_leaf_set_sizes_unlabeled_symmetric_difference(other);
}

unsigned long DistanceTree::get_unweighted_labeled_symmetric_difference(DistanceTree & other) {
    throw std::logic_error("Not Implemented: DistanceTree::get_unweighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long DistanceTree::get_weighted_labeled_symmetric_difference(DistanceTree & other) {
    throw std::logic_error("Not Implemented: DistanceTree::get_weighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long DistanceTree::calc_set_symmetric_difference(
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
} // namespace pstrudel
