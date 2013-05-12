#include "distancetree.hpp"

namespace pstrudel {

//////////////////////////////////////////////////////////////////////////////
// PairwiseTipDistanceProfileCalculator

PairwiseTipDistanceProfileCalculator & PairwiseTipDistanceProfileCalculator::operator=(const PairwiseTipDistanceProfileCalculator & other) {
    this->unweighted_pairwise_tip_distance_profile_ = other.unweighted_pairwise_tip_distance_profile_;
    this->weighted_pairwise_tip_distance_profile_ = other.weighted_pairwise_tip_distance_profile_;
    return *this;
}

void PairwiseTipDistanceProfileCalculator::build_pairwise_tip_distance_profiles() {
    std::multiset<unsigned long> unweighted_distances;
    std::multiset<double> weighted_distances;
    std::vector<double> scaled_weighted_distances;
    auto f = [&unweighted_distances, &weighted_distances] (DistanceNodeValue &,
            DistanceNodeValue &,
            unsigned long u,
            double w) {
        unweighted_distances.insert(u);
        weighted_distances.insert(w);
    };
    this->tree_.calc_pairwise_tip_distances(f);
    this->unweighted_pairwise_tip_distance_profile_.set_data(unweighted_distances.begin(),
            unweighted_distances.end());
    this->weighted_pairwise_tip_distance_profile_.set_data(weighted_distances.begin(),
            weighted_distances.end());
    double tree_length = this->tree_.get_total_tree_length();
    if (tree_length > 0) {
        scaled_weighted_distances.reserve(weighted_distances.size());
        for (auto & wd : weighted_distances) {
            scaled_weighted_distances.push_back(wd / tree_length);
        }
    } else {
        scaled_weighted_distances.insert(scaled_weighted_distances.end(), weighted_distances.size(), 0.0);
    }
    this->scaled_weighted_pairwise_tip_distance_profile_.set_data(scaled_weighted_distances.begin(),
            scaled_weighted_distances.end());
}

double PairwiseTipDistanceProfileCalculator::get_unweighted_distance(PairwiseTipDistanceProfileCalculator & other) {
    if (this->unweighted_pairwise_tip_distance_profile_.empty()) {
        this->build_pairwise_tip_distance_profiles();
    }
    if (other.unweighted_pairwise_tip_distance_profile_.empty()) {
        other.build_pairwise_tip_distance_profiles();
    }
    return this->unweighted_pairwise_tip_distance_profile_.get_distance(other.unweighted_pairwise_tip_distance_profile_);
}

double PairwiseTipDistanceProfileCalculator::get_weighted_distance(PairwiseTipDistanceProfileCalculator & other) {
    if (this->weighted_pairwise_tip_distance_profile_.empty()) {
        this->build_pairwise_tip_distance_profiles();
    }
    if (other.weighted_pairwise_tip_distance_profile_.empty()) {
        other.build_pairwise_tip_distance_profiles();
    }
    return this->weighted_pairwise_tip_distance_profile_.get_distance(other.weighted_pairwise_tip_distance_profile_);
}

double PairwiseTipDistanceProfileCalculator::get_scaled_weighted_distance(PairwiseTipDistanceProfileCalculator & other) {
    if (this->scaled_weighted_pairwise_tip_distance_profile_.empty()) {
        this->build_pairwise_tip_distance_profiles();
    }
    if (other.scaled_weighted_pairwise_tip_distance_profile_.empty()) {
        other.build_pairwise_tip_distance_profiles();
    }
    return this->scaled_weighted_pairwise_tip_distance_profile_.get_distance(other.scaled_weighted_pairwise_tip_distance_profile_);
}

//////////////////////////////////////////////////////////////////////////////
// SymmetricDifferenceCalculator

SymmetricDifferenceCalculator & SymmetricDifferenceCalculator::operator=(const SymmetricDifferenceCalculator & other) {
    this->subtree_leaf_set_sizes_ = other.subtree_leaf_set_sizes_;
    return *this;
}

void SymmetricDifferenceCalculator::calc_subtree_leaf_set_sizes() {
    this->subtree_leaf_set_sizes_.clear();
    for (auto ndi = this->tree_.postorder_begin(); ndi != this->tree_.postorder_end(); ++ndi) {
        if (ndi.is_leaf()) {
            ndi->set_num_leaves(1);
        } else {
            unsigned long nleaves = 0;
            unsigned long ndesc = 0;
            for (auto chi = this->tree_.children_begin(ndi); chi != this->tree_.children_end(ndi); ++chi) {
                nleaves += chi->get_num_leaves();
                ndesc += 1;
            }
            ndi->set_num_leaves(nleaves);
            this->subtree_leaf_set_sizes_.insert(nleaves);
        }
    }
}

unsigned long SymmetricDifferenceCalculator::calc_leaf_set_sizes_unlabeled_symmetric_difference(SymmetricDifferenceCalculator & other) {
    if (this->subtree_leaf_set_sizes_.empty()) {
        this->calc_subtree_leaf_set_sizes();
    }
    if (other.subtree_leaf_set_sizes_.empty()) {
        other.calc_subtree_leaf_set_sizes();
    }
    return SymmetricDifferenceCalculator::calc_set_symmetric_difference(this->subtree_leaf_set_sizes_, other.subtree_leaf_set_sizes_);
}

unsigned long SymmetricDifferenceCalculator::get_unlabeled_symmetric_difference(SymmetricDifferenceCalculator & other) {
    return this->calc_leaf_set_sizes_unlabeled_symmetric_difference(other);
}

unsigned long SymmetricDifferenceCalculator::get_unweighted_labeled_symmetric_difference(SymmetricDifferenceCalculator & other) {
    throw std::logic_error("Not Implemented: SymmetricDifferenceCalculator::get_unweighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long SymmetricDifferenceCalculator::get_weighted_labeled_symmetric_difference(SymmetricDifferenceCalculator & other) {
    throw std::logic_error("Not Implemented: SymmetricDifferenceCalculator::get_weighted_labeled_symmetric_difference()");
    return 0;
}

unsigned long SymmetricDifferenceCalculator::calc_set_symmetric_difference(
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

//////////////////////////////////////////////////////////////////////////////
// DistanceTree

DistanceTree::DistanceTree(bool is_rooted)
        : platypus::StandardTree<DistanceNodeValue>(is_rooted)
          , number_of_tips_(0)
          , total_tree_length_(0.0)
          , pairwise_tip_distance_profile_calculator_(*this)
          , symmetric_difference_calculator_(*this) {
}

DistanceTree::DistanceTree(DistanceTree && other)
        : platypus::StandardTree<DistanceNodeValue>(other)
        , number_of_tips_(other.number_of_tips_)
        , total_tree_length_(other.total_tree_length_)
        , pairwise_tip_distance_profile_calculator_(*this)
        , symmetric_difference_calculator_(*this) {
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
}

DistanceTree::DistanceTree(const DistanceTree & other)
        : platypus::StandardTree<DistanceNodeValue>(other)
        , number_of_tips_(other.number_of_tips_)
        , total_tree_length_(other.total_tree_length_)
        , pairwise_tip_distance_profile_calculator_(*this)
        , symmetric_difference_calculator_(*this) {
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
}

DistanceTree & DistanceTree::operator=(const DistanceTree & other) {
    platypus::StandardTree<DistanceNodeValue>::operator=(other);
    this->number_of_tips_ = other.number_of_tips_;
    this->total_tree_length_ = other.total_tree_length_;
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
    return *this;
}

} // namespace pstrudel
