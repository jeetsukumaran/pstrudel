#ifndef PSTRUDEL_DISTANCE_TREE_HPP
#define PSTRUDEL_DISTANCE_TREE_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <colugo/utility.hpp>
#include <platypus/model/standardtree.hpp>
#include "profile.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// DistanceNodeValue

class DistanceNodeValue : public platypus::StandardNodeValue {
    public:
        DistanceNodeValue()
            : num_leaves_(0) {
        }
        DistanceNodeValue(DistanceNodeValue&& other)
            : platypus::StandardNodeValue(std::move(other))
            , desc_path_lens_(std::move(other.desc_path_lens_))
            , num_leaves_(other.num_leaves_) {
        }
        DistanceNodeValue(const DistanceNodeValue& other)
            : platypus::StandardNodeValue(other)
            , desc_path_lens_(other.desc_path_lens_)
            , num_leaves_(other.num_leaves_) {
        }
        DistanceNodeValue& operator=(const DistanceNodeValue& other) {
            platypus::StandardNodeValue::operator=(other);
            this->desc_path_lens_ = other.desc_path_lens_;
            this->num_leaves_ = other.num_leaves_;
            return *this;
        }
        void clear() override {
            platypus::StandardNodeValue::clear();
            this->desc_path_lens_.clear();
            this->num_leaves_ = 0;
        }
        void set_desc_path_len(DistanceNodeValue & nd, double len, unsigned long step_count) {
            this->desc_path_lens_[&nd] = std::make_pair(len, step_count);
        }
        void set_desc_path_len(DistanceNodeValue * nd, double len, unsigned long step_count) {
            this->desc_path_lens_[nd] = std::make_pair(len, step_count);
        }
        double get_desc_path_weight(DistanceNodeValue * nd) {
            return this->desc_path_lens_[nd].first;
        }
        unsigned long get_desc_path_count(DistanceNodeValue * nd) {
            return this->desc_path_lens_[nd].second;
        }
        std::map<DistanceNodeValue *, std::pair<double, unsigned long>>& get_desc_paths() {
            return this->desc_path_lens_;
        }
        unsigned long get_num_leaves() const {
            return this->num_leaves_;
        }
        void set_num_leaves(unsigned long num_leaves) {
            this->num_leaves_ = num_leaves;
        }
    private:
        std::map<DistanceNodeValue *, std::pair<double, unsigned long>>   desc_path_lens_;
        unsigned long                                                     num_leaves_;
}; // DistanceNodeValue

////////////////////////////////////////////////////////////////////////////////
// DistanceTree

class DistanceTree : public platypus::StandardTree<DistanceNodeValue> {

    public:
        typedef std::unordered_multiset<unsigned long> SizesSetType;
        const char * UNWEIGHTED_PAIRWISE_TIP = "unweighted pairwise tip";
        const char * WEIGHTED_PAIRWISE_TIP = "weighted pairwise tip";

    public:
        DistanceTree(bool is_rooted=true);
        DistanceTree(DistanceTree && other);
        DistanceTree(const DistanceTree & other);
        DistanceTree & operator=(const DistanceTree & other);
        ~DistanceTree() override {}

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

        double get_weighted_pairwise_tip_distance(DistanceNodeValue& tip1, DistanceNodeValue& tip2) {
            if (&tip1 == &tip2) {
                return 0.0;
            }
            auto v = this->find_pairwise_tip_distance(this->weighted_pairwise_tip_distance_, tip1, tip2);
            return v->second;
        }

        unsigned long get_unweighted_pairwise_tip_distance(DistanceNodeValue& tip1, DistanceNodeValue& tip2) {
            if (&tip1 == &tip2) {
                return 0;
            }
            auto v = this->find_pairwise_tip_distance(this->unweighted_pairwise_tip_distance_, tip1, tip2);
            return v->second;
        }

        double get_unweighted_subprofile_distance(const DistanceTree& other) const {
            return this->profile_.calc_subprofile_distance(other.profile_, this->UNWEIGHTED_PAIRWISE_TIP);
        }

        double get_weighted_subprofile_distance(const DistanceTree& other) const {
            return this->profile_.calc_subprofile_distance(other.profile_, this->WEIGHTED_PAIRWISE_TIP);
        }

        void calc_subtree_sizes();
        unsigned long calc_leaf_set_sizes_unlabeled_symmetric_difference(DistanceTree & other);
        unsigned long get_unlabeled_symmetric_difference(DistanceTree & other);
        unsigned long get_unweighted_labeled_symmetric_difference(DistanceTree & other);
        unsigned long get_weighted_labeled_symmetric_difference(DistanceTree & other);

    public:

        static unsigned long calc_set_symmetric_difference(
                const SizesSetType & set1,
                const SizesSetType & set2,
                SizesSetType * common = nullptr,
                SizesSetType * unmatched1 = nullptr,
                SizesSetType * unmatched2 = nullptr);

    private:
        template <typename T>
        typename T::const_iterator find_pairwise_tip_distance(T& container, DistanceNodeValue& tip1, DistanceNodeValue& tip2) {
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
        std::map< std::pair<DistanceNodeValue *, DistanceNodeValue *>, double >             weighted_pairwise_tip_distance_;
        std::map< std::pair<DistanceNodeValue *, DistanceNodeValue *>, unsigned long >      unweighted_pairwise_tip_distance_;
        unsigned long                                                                       number_of_tips_;
        double                                                                              total_tree_length_;
        Profile                                                                             profile_;
        SizesSetType                                                                        subtree_leaf_set_sizes_;

}; // DistanceTree

} // namespace pstrudel
#endif
