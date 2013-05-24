#include <map>
#include <platypus/model/coalescent.hpp>
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
    std::vector<unsigned long> unweighted_distances;
    std::vector<double> weighted_distances;
    std::vector<double> scaled_weighted_distances;
    unsigned long nc = (this->tree_.get_num_tips() * (this->tree_.get_num_tips() - 1))/2;
    weighted_distances.reserve(nc);
    unweighted_distances.reserve(nc);
    scaled_weighted_distances.reserve(nc);
    auto f = [&unweighted_distances, &weighted_distances] (DistanceNodeValue &,
            DistanceNodeValue &,
            unsigned long u,
            double w) {
        unweighted_distances.push_back(u);
        weighted_distances.push_back(w);
    };
    this->tree_.calc_pairwise_tip_distances(f);
    this->unweighted_pairwise_tip_distance_profile_.set_data(unweighted_distances.begin(),
            unweighted_distances.end());
    this->weighted_pairwise_tip_distance_profile_.set_data(weighted_distances.begin(),
            weighted_distances.end());
    double tree_length = this->tree_.get_total_tree_length();
    if (tree_length > 0) {
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
// LineageThroughTimeProfileCalculator

LineageThroughTimeProfileCalculator & LineageThroughTimeProfileCalculator::operator=(const LineageThroughTimeProfileCalculator & other) {
    this->lineage_accumulation_through_time_profile_ = other.lineage_accumulation_through_time_profile_;
    this->lineage_splitting_time_profile_ = other.lineage_splitting_time_profile_;
    this->scaled_lineage_splitting_time_profile_ = other.scaled_lineage_splitting_time_profile_;
    this->max_leaf_distance_ = other.max_leaf_distance_;
    return *this;
}

double LineageThroughTimeProfileCalculator::get_lineage_accumulation_profile_distance(LineageThroughTimeProfileCalculator & other) {
    if (this->lineage_accumulation_through_time_profile_.empty()) {
        this->build_lineage_accumulation_through_time_profile();
    }
    if (other.lineage_accumulation_through_time_profile_.empty()) {
        other.build_lineage_accumulation_through_time_profile();
    }
    return this->lineage_accumulation_through_time_profile_.get_distance(other.lineage_accumulation_through_time_profile_);
}

double LineageThroughTimeProfileCalculator::get_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other) {
    if (this->lineage_splitting_time_profile_.empty()) {
        this->build_lineage_splitting_time_profile();
    }
    if (other.lineage_splitting_time_profile_.empty()) {
        other.build_lineage_splitting_time_profile();
    }
    return this->lineage_splitting_time_profile_.get_distance(other.lineage_splitting_time_profile_);
}

double LineageThroughTimeProfileCalculator::get_scaled_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other) {
    if (this->scaled_lineage_splitting_time_profile_.empty()) {
        this->build_lineage_splitting_time_profile();
    }
    if (other.scaled_lineage_splitting_time_profile_.empty()) {
        other.build_lineage_splitting_time_profile();
    }
    return this->scaled_lineage_splitting_time_profile_.get_distance(other.scaled_lineage_splitting_time_profile_);
}

void LineageThroughTimeProfileCalculator::clear() {
    this->max_leaf_distance_ = 0.0;
    this->lineage_accumulation_through_time_profile_.clear();
    this->lineage_splitting_time_profile_.clear();
    this->scaled_lineage_splitting_time_profile_.clear();
}

unsigned long LineageThroughTimeProfileCalculator::get_default_num_transects() {
    if (this->tree_.get_num_tips() == 0) {
        this->tree_.calc_num_tips();
    }
    return (this->tree_.get_num_tips() - 1) * 10;
}

void LineageThroughTimeProfileCalculator::calc_node_root_distances() {
    this->max_leaf_distance_ = 0.0;
    double distance = 0.0;
    for (auto ndi = this->tree_.preorder_begin(); ndi != this->tree_.preorder_end(); ++ndi) {
        if (ndi.parent_node() == nullptr) {
            distance = 0.0;
        } else {
            distance = ndi.parent().get_root_distance() + ndi->get_edge_length();
        }
        ndi->set_root_distance(distance);
        if (distance > this->max_leaf_distance_) {
            this->max_leaf_distance_ = distance;
        }
    }
}

std::vector<double> LineageThroughTimeProfileCalculator::build_transect_offsets(unsigned long num_transects) {
    std::vector<double> transect_offsets;
    if (this->max_leaf_distance_ == 0) {
        this->calc_node_root_distances();
    }
    double transect_step = (static_cast<double>(this->max_leaf_distance_) / num_transects);
    // for (double t = transect_step; t < this->max_leaf_distance_; t += transect_step) {
    double offset = transect_step;
    for (unsigned int i; i <= num_transects; ++i) {
        if (offset > this->max_leaf_distance_) {
            transect_offsets.push_back(this->max_leaf_distance_);
            break;
        }
        transect_offsets.push_back(offset);
        offset += transect_step;
    }
    return transect_offsets;
}

const Profile & LineageThroughTimeProfileCalculator::build_lineage_accumulation_through_time_profile(const std::vector<double> & transect_offsets) {
    this->lineage_accumulation_through_time_profile_.clear();
    if (this->tree_.get_num_tips() == 0) {
        this->tree_.calc_num_tips();
    }
    double num_tips = this->tree_.get_num_tips();
    if (this->max_leaf_distance_ == 0) {
        this->calc_node_root_distances();
    }
    this->calc_node_root_distances();
    std::vector<double> num_lineages(transect_offsets.size(), 0.0);
    for (auto ndi = this->tree_.preorder_begin(); ndi != this->tree_.preorder_end(); ++ndi) {
        if (ndi.parent_node() != nullptr) {
            auto transect_offsets_begin = transect_offsets.begin();
            for (auto toi = transect_offsets_begin ; toi != transect_offsets.end(); ++toi) {
                if (ndi->get_root_distance() >= *toi && ndi.parent().get_root_distance() < *toi) {
                    auto idx = toi - transect_offsets_begin;
                    num_lineages[idx] += 1.0/num_tips;
                }
            }
        }
    }
    this->lineage_accumulation_through_time_profile_.set_data(num_lineages.begin(), num_lineages.end(), false);
    return this->lineage_accumulation_through_time_profile_;
}

const Profile & LineageThroughTimeProfileCalculator::build_lineage_accumulation_through_time_profile(unsigned long num_transects) {
    if (num_transects == 0) {
        num_transects = this->get_default_num_transects();
    }
    auto transect_offsets = this->build_transect_offsets(num_transects);
    return this->build_lineage_accumulation_through_time_profile(transect_offsets);
}

std::pair<const Profile &, const Profile &> LineageThroughTimeProfileCalculator::build_lineage_splitting_time_profile() {
    std::pair<const Profile &, const Profile &> ret_val(this->lineage_splitting_time_profile_, this->scaled_lineage_splitting_time_profile_);
    this->lineage_splitting_time_profile_.clear();
    this->scaled_lineage_splitting_time_profile_.clear();
    if (this->tree_.get_num_tips() == 0) {
        this->tree_.calc_num_tips();
    }
    double num_tips = this->tree_.get_num_tips();
    if (this->max_leaf_distance_ == 0) {
        this->calc_node_root_distances();
    }
    if (this->tree_.get_total_tree_length() == 0) {
        this->tree_.calc_total_tree_length();
    }
    std::vector<double> splitting_times;
    std::vector<double> scaled_splitting_times;
    splitting_times.reserve(num_tips);
    scaled_splitting_times.reserve(num_tips);
    double tree_length = this->tree_.get_total_tree_length();
    if (tree_length > 0.0) {
        double dist = 0.0;
        for (auto ndi = this->tree_.preorder_begin(); ndi != this->tree_.preorder_end(); ++ndi) {
            if (!ndi.is_leaf()) {
                dist = ndi->get_root_distance();
                splitting_times.push_back(dist);
                scaled_splitting_times.push_back(dist/tree_length);
            }
        }
    } else {
        splitting_times.insert(splitting_times.end(), num_tips-1, 0);
        scaled_splitting_times.insert(scaled_splitting_times.end(), num_tips-1, 0);
    }
    this->lineage_splitting_time_profile_.set_data(splitting_times.begin(), splitting_times.end(), true);
    this->scaled_lineage_splitting_time_profile_.set_data(scaled_splitting_times.begin(), scaled_splitting_times.end(), true);
    return ret_val;
}

//////////////////////////////////////////////////////////////////////////////
// DistanceTree

const std::vector<std::string> DistanceTree::tree_pattern_y_distance_names_{
        "y.ptd.uw",
        "y.ptd.wt",
        "y.ltt",
        "y.lst",
}; // static cons ttree_pattern_names_

DistanceTree::DistanceTree(bool is_rooted)
        : platypus::StandardTree<DistanceNodeValue>(is_rooted)
          , number_of_tips_(0)
          , total_tree_length_(0.0)
          , pairwise_tip_distance_profile_calculator_(*this)
          , symmetric_difference_calculator_(*this)
          , lineage_through_time_calculator_(*this) {
}

DistanceTree::DistanceTree(DistanceTree && other)
        : platypus::StandardTree<DistanceNodeValue>(other)
         , number_of_tips_(other.number_of_tips_)
         , total_tree_length_(other.total_tree_length_)
         , pairwise_tip_distance_profile_calculator_(*this)
         , symmetric_difference_calculator_(*this)
         , lineage_through_time_calculator_(*this) {
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
    this->lineage_through_time_calculator_ = other.lineage_through_time_calculator_;
}

DistanceTree::DistanceTree(const DistanceTree & other)
        : platypus::StandardTree<DistanceNodeValue>(other)
         , number_of_tips_(other.number_of_tips_)
         , total_tree_length_(other.total_tree_length_)
         , pairwise_tip_distance_profile_calculator_(*this)
         , symmetric_difference_calculator_(*this)
         , lineage_through_time_calculator_(*this) {
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
    this->lineage_through_time_calculator_ = other.lineage_through_time_calculator_;
}

DistanceTree & DistanceTree::operator=(const DistanceTree & other) {
    platypus::StandardTree<DistanceNodeValue>::operator=(other);
    this->number_of_tips_ = other.number_of_tips_;
    this->total_tree_length_ = other.total_tree_length_;
    this->pairwise_tip_distance_profile_calculator_ = other.pairwise_tip_distance_profile_calculator_;
    this->symmetric_difference_calculator_ = other.symmetric_difference_calculator_;
    this->lineage_through_time_calculator_ = other.lineage_through_time_calculator_;
    return *this;
}

std::vector<DistanceTree::node_type *> DistanceTree::get_nodes_in_level_order() {
    std::vector<DistanceTree::node_type *> nodes_in_level_order;
    std::map<DistanceTree::node_type *, unsigned long> steps_from_root;
    for (auto ndi = this->preorder_begin(); ndi != this->preorder_end(); ++ndi) {
        DistanceTree::node_type * nptr = ndi.node();
        if (ndi.parent_node() == nullptr) {
            steps_from_root[nptr] = 0;
        } else {
            steps_from_root[nptr] = steps_from_root[ndi.parent_node()] + 1;
        }
    }
    std::multimap<unsigned long, DistanceTree::node_type *> steps_from_root_inv;
    for (auto & si : steps_from_root) {
        steps_from_root_inv.insert(std::make_pair(si.second, si.first));
    }
    nodes_in_level_order.reserve(steps_from_root_inv.size());
    for (auto & si : steps_from_root_inv) {
        nodes_in_level_order.push_back(si.second);
    }
    return nodes_in_level_order;
}

// type: 0 = mean coalescent, 1 = random, 2 = uniform 3 = anti-coalescent
void DistanceTree::add_edge_lengths(int regime) {
    auto nodes_in_level_order = this->get_nodes_in_level_order();
    if (this->number_of_tips_ == 0) {
        for (auto ndi = this->leaf_begin(); ndi != this->leaf_end(); ++ndi, ++this->number_of_tips_) {
        }
    }
    std::map<DistanceTree::node_type *, double> node_ages;
    unsigned long num_lineages = this->number_of_tips_;
    unsigned long total_lineages = num_lineages;
    double prev_age = 0.0;
    for (auto ndi = nodes_in_level_order.rbegin(); ndi != nodes_in_level_order.rend(); ++ndi) {
        if ((*ndi)->is_leaf()) {
            node_ages[*ndi] = 0.0;
        } else {
            double wt;
            if (regime == 0) {
                wt = platypus::coalescent::expected_time_to_coalescence(num_lineages, 1.0, 2);
                // std::cerr << num_lineages << ", " << wt << std::endl;
            } else if (regime == 1) {
                wt = platypus::coalescent::expected_time_to_coalescence(total_lineages - num_lineages + 2, 1.0, 2);
            } else if (regime == 2) {
                wt = 1.0;
            } else {
                throw std::runtime_error("Unsupported regime: " + std::to_string(regime));
            }
            node_ages[*ndi] = prev_age + wt;
            prev_age = node_ages[*ndi];
            --num_lineages;
        }
    }
    for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
        if (ndi.parent_node() != nullptr) {
            ndi->set_edge_length(node_ages[ndi.parent_node()] - node_ages[ndi.node()]);
        } else {
            ndi->set_edge_length(0.0);
        }
    }
}

} // namespace pstrudel
