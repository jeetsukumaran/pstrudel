#include <algorithm>
#include <iostream>
#include "pairwise_distance_tree.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// PairwiseDistanceTree

PairwiseDistanceTree::PairwiseDistanceTree(bool is_rooted)
        : BasicTree(is_rooted) {
}

PairwiseDistanceTree::PairwiseDistanceTree(PairwiseDistanceTree && other)
        : BasicTree(other)
          , weighted_pairwise_tip_distance_(std::move(other.weighted_pairwise_tip_distance_))
          , unweighted_pairwise_tip_distance_(std::move(other.unweighted_pairwise_tip_distance_))
          , total_tree_length_(other.total_tree_length_)
          , number_of_tips_(other.number_of_tips_)
          , profile_(other.profile_) {
}

PairwiseDistanceTree::PairwiseDistanceTree(const PairwiseDistanceTree & other)
        : BasicTree(other)
          , weighted_pairwise_tip_distance_(other.weighted_pairwise_tip_distance_)
          , unweighted_pairwise_tip_distance_(other.unweighted_pairwise_tip_distance_)
          , total_tree_length_(other.total_tree_length_)
          , number_of_tips_(other.number_of_tips_)
          , profile_(other.profile_) {
}

PairwiseDistanceTree & PairwiseDistanceTree::operator=(const PairwiseDistanceTree & other) {
    BasicTree::operator=(other);
    this->weighted_pairwise_tip_distance_ = other.weighted_pairwise_tip_distance_;
    this->unweighted_pairwise_tip_distance_ = other.unweighted_pairwise_tip_distance_;
    this->total_tree_length_ = other.total_tree_length_;
    this->number_of_tips_ = other.number_of_tips_;
    this->profile_ = other.profile_;
    return *this;
}

void PairwiseDistanceTree::calc_profile_metrics() {
    this->profile_.clear();
    this->calc_pairwise_tip_distances();
}

void PairwiseDistanceTree::build_profile() {
    this->profile_.build_subprofiles();
}

void PairwiseDistanceTree::calc_pairwise_tip_distances() {
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

void PairwiseDistanceTree::poll_max_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    this->profile_.poll_max_num_interpolated_points(num_interpolated_points);
}

void PairwiseDistanceTree::set_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    this->profile_.set_num_interpolated_points(num_interpolated_points);
}

} // pstrudel
