#ifndef PSTRUDEL_PAIRWISE_DISTANCE_TREE_HPP
#define PSTRUDEL_PAIRWISE_DISTANCE_TREE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include "colugo-utilities/src/utility.hpp"
#include "platypus-phyloinformary/src/tree.hpp"
#include "basic_tree.hpp"
#include "profile.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// PairwiseDistanceNodeValue

class PairwiseDistanceNodeValue : public BasicNodeValue {
    public:
        PairwiseDistanceNodeValue(int index = -1)
            : BasicNodeValue(0.0) {
        }
        PairwiseDistanceNodeValue(PairwiseDistanceNodeValue&& other)
            : BasicNodeValue(std::move(other))
              , desc_path_lens_(std::move(other.desc_path_lens_)) {
        }
        PairwiseDistanceNodeValue(const PairwiseDistanceNodeValue& other)
            : BasicNodeValue(other)
              , desc_path_lens_(other.desc_path_lens_) {
        }
        PairwiseDistanceNodeValue& operator=(const PairwiseDistanceNodeValue& other) {
            BasicNodeValue::operator=(other);
            this->desc_path_lens_ = other.desc_path_lens_;
            return *this;
        }
        void clear() override {
            BasicNodeValue::clear();
            this->desc_path_lens_.clear();
        }
        void set_desc_path_len(PairwiseDistanceNodeValue & nd, double len, unsigned long step_count) {
            this->desc_path_lens_[&nd] = std::make_pair(len, step_count);
        }
        void set_desc_path_len(PairwiseDistanceNodeValue * nd, double len, unsigned long step_count) {
            this->desc_path_lens_[nd] = std::make_pair(len, step_count);
        }
        double get_desc_path_weight(PairwiseDistanceNodeValue * nd) {
            return this->desc_path_lens_[nd].first;
        }
        unsigned long get_desc_path_count(PairwiseDistanceNodeValue * nd) {
            return this->desc_path_lens_[nd].second;
        }
        std::map<PairwiseDistanceNodeValue *, std::pair<double, unsigned long>>& get_desc_paths() {
            return this->desc_path_lens_;
        }
    private:
        int             index_;
        double          edge_length_;
        std::string     label_;
        std::map<PairwiseDistanceNodeValue *, std::pair<double, unsigned long>>    desc_path_lens_;
}; // PairwiseDistanceNodeValue

////////////////////////////////////////////////////////////////////////////////
// PairwiseDistanceTree

class PairwiseDistanceTree : public BasicTree<PairwiseDistanceNodeValue> {

    public:
        const char * UNWEIGHTED_PAIRWISE_TIP = "unweighted pairwise tip";
        const char * WEIGHTED_PAIRWISE_TIP = "weighted pairwise tip";

    public:
        PairwiseDistanceTree(bool is_rooted=true);
        PairwiseDistanceTree(PairwiseDistanceTree && other);
        PairwiseDistanceTree(const PairwiseDistanceTree & other);
        PairwiseDistanceTree & operator=(const PairwiseDistanceTree & other);
        ~PairwiseDistanceTree() override {}

        /////////////////////////////////////////////////////////////////////////
        // profile management

        void clear();
        void calc_profile_metrics();
        void build_profile();
        void calc_pairwise_tip_distances();
        void poll_max_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points);
        void set_num_interpolated_profile_points(std::map<const std::string, unsigned long>& num_interpolated_points);

        /////////////////////////////////////////////////////////////////////////
        // querying metrics

        double get_weighted_pairwise_tip_distance(PairwiseDistanceNodeValue& tip1, PairwiseDistanceNodeValue& tip2) {
            if (&tip1 == &tip2) {
                return 0.0;
            }
            auto v = this->find_pairwise_tip_distance(this->weighted_pairwise_tip_distance_, tip1, tip2);
            return v->second;
        }

        unsigned long get_unweighted_pairwise_tip_distance(PairwiseDistanceNodeValue& tip1, PairwiseDistanceNodeValue& tip2) {
            if (&tip1 == &tip2) {
                return 0;
            }
            auto v = this->find_pairwise_tip_distance(this->unweighted_pairwise_tip_distance_, tip1, tip2);
            return v->second;
        }

        double get_unweighted_subprofile_distance(const PairwiseDistanceTree& other) const {
            return this->profile_.calc_subprofile_distance(other.profile_, this->UNWEIGHTED_PAIRWISE_TIP);
        }

        double get_weighted_subprofile_distance(const PairwiseDistanceTree& other) const {
            return this->profile_.calc_subprofile_distance(other.profile_, this->WEIGHTED_PAIRWISE_TIP);
        }

        double get_distance(const PairwiseDistanceTree& other) const {
            return this->profile_.calc_distance(other.profile_);
        }

    private:
        template <typename T>
        typename T::const_iterator find_pairwise_tip_distance(T& container, PairwiseDistanceNodeValue& tip1, PairwiseDistanceNodeValue& tip2) {
            auto v = container.find(std::make_pair(&tip1, &tip2));
            if (v == container.end()) {
                v = container.find(std::make_pair(&tip2, &tip1));
                if (v == container.end()) {
                    colugo::colugo_abort("Nodes not found on tree");
                }
                return v;
            }
            return v;
        }

    private:
        std::map< std::pair<PairwiseDistanceNodeValue *, PairwiseDistanceNodeValue *>, double >           weighted_pairwise_tip_distance_;
        std::map< std::pair<PairwiseDistanceNodeValue *, PairwiseDistanceNodeValue *>, unsigned long >    unweighted_pairwise_tip_distance_;
        double                                                                          total_tree_length_;
        unsigned long                                                                   number_of_tips_;
        Profile                                                                         profile_;

}; // PairwiseDistanceTree


} // namespace pstrudel

#endif
