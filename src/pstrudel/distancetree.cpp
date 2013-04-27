#include "distancetree.hpp"

namespace pstrudel {

DistanceTree::DistanceTree(bool is_rooted)
        : platypus::StandardTree<DistanceNodeValue>(is_rooted) {
}

DistanceTree::DistanceTree(DistanceTree && other)
    : platypus::StandardTree<DistanceNodeValue>(other)
    , weighted_pairwise_tip_distance_(std::move(other.weighted_pairwise_tip_distance_))
    , unweighted_pairwise_tip_distance_(std::move(other.unweighted_pairwise_tip_distance_))
    , number_of_tips_(other.number_of_tips_)
    , total_tree_length_(other.total_tree_length_)
    , profile_(std::move(other.profile_))
    , subtree_leaf_set_sizes_(std::move(other.subtree_leaf_set_sizes_)) {
}

DistanceTree::DistanceTree(const DistanceTree & other)
    : platypus::StandardTree<DistanceNodeValue>(other)
    , weighted_pairwise_tip_distance_(other.weighted_pairwise_tip_distance_)
    , unweighted_pairwise_tip_distance_(other.unweighted_pairwise_tip_distance_)
    , number_of_tips_(other.number_of_tips_)
    , total_tree_length_(other.total_tree_length_)
    , profile_(other.profile_)
    , subtree_leaf_set_sizes_(other.subtree_leaf_set_sizes_) {
}

DistanceTree & DistanceTree::operator=(const DistanceTree & other) {
    platypus::StandardTree<DistanceNodeValue>::operator=(other);
    this->weighted_pairwise_tip_distance_ = other.weighted_pairwise_tip_distance_;
    this->unweighted_pairwise_tip_distance_ = other.unweighted_pairwise_tip_distance_;
    this->number_of_tips_ = other.number_of_tips_;
    this->total_tree_length_ = other.total_tree_length_;
    this->profile_ = other.profile_;
    this->subtree_leaf_set_sizes_ = other.subtree_leaf_set_sizes_;
    return *this;
}

void DistanceTree::calc_profile_metrics() {
    this->profile_.clear();
    this->calc_pairwise_tip_distances();
}

void DistanceTree::build_profile() {
    this->profile_.build_subprofiles();
}

void DistanceTree::calc_pairwise_tip_distances() {
    this->weighted_pairwise_tip_distance_.clear();
    this->unweighted_pairwise_tip_distance_.clear();
    this->total_tree_length_ = 0.0;
    this->number_of_tips_ = 0;
    std::multiset<double> weighted_distances;
    std::multiset<unsigned long> unweighted_distances;
    for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        this->total_tree_length_ += ndi->get_edge_length();
        if (ndi.is_leaf()) {
            ++this->number_of_tips_;
            ndi->set_desc_path_len(*ndi, 0.0, 0);
        } else {
            for (auto chi1 = this->children_begin(ndi) ; chi1 != this->children_end(ndi) ; ++chi1) {
                for (auto & desc1 : chi1->get_desc_paths()) {
                    auto & desc1_nd = desc1.first;
                    auto & desc1_weight = desc1.second.first;
                    auto & desc1_steps = desc1.second.second;
                    ndi->set_desc_path_len(desc1_nd, chi1->get_edge_length() + desc1_weight, desc1_steps + 1);
                    for (auto chi2 = this->next_sibling_begin(chi1) ; chi2 != this->next_sibling_end(chi1) ; ++chi2) {
                        for (auto & desc2 : chi2->get_desc_paths()) {
                            auto & desc2_nd = desc2.first;
                            auto & desc2_weight = desc2.second.first;
                            auto & desc2_steps = desc2.second.second;
                            auto weight = ndi->get_desc_path_weight(desc1_nd) + desc2_weight + chi2->get_edge_length();
                            auto steps = ndi->get_desc_path_count(desc1_nd) + desc2_steps + 1;
                            this->weighted_pairwise_tip_distance_[std::make_pair(&(*desc1_nd), &(*desc2_nd))] = weight;
                            this->unweighted_pairwise_tip_distance_[std::make_pair(&(*desc1_nd), &(*desc2_nd))] = steps;
                            weighted_distances.insert(weight);
                            unweighted_distances.insert(steps);
                        }
                    }
                }
            }
        }
    }
    this->profile_.add_subprofile_source(this->UNWEIGHTED_PAIRWISE_TIP,
            unweighted_distances.begin(),
            unweighted_distances.end(),
            unweighted_distances.size(),
            false);
    this->profile_.add_subprofile_source(this->WEIGHTED_PAIRWISE_TIP,
            weighted_distances.begin(),
            weighted_distances.end(),
            weighted_distances.size(),
            false);
}

void DistanceTree::poll_max_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    this->profile_.poll_max_num_interpolated_points(num_interpolated_points);
}

void DistanceTree::set_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    this->profile_.set_num_interpolated_points(num_interpolated_points);
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
