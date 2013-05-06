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

        /////////////////////////////////////////////////////////////////////////
        // lifecycle and assignment

        DistanceTree(bool is_rooted=true);
        DistanceTree(DistanceTree && other);
        DistanceTree(const DistanceTree & other);
        DistanceTree & operator=(const DistanceTree & other);
        ~DistanceTree() override {}

        /////////////////////////////////////////////////////////////////////////
        // basic stats

        void set_num_tips(unsigned long n) {
            this->number_of_tips_ = n;
        }
        unsigned long get_num_tips() const {
            return this->number_of_tips_;
        }
        void set_total_tree_length(double v) {
            this->total_tree_length_ = v;
        }
        double get_total_tree_length() const {
            return this->total_tree_length_;
        }

        /////////////////////////////////////////////////////////////////////////
        // profile management

        void build_pairwise_tip_distance_profiles() {
            std::multiset<unsigned long> unweighted_distances;
            std::multiset<double> weighted_distances;
            auto f = [&unweighted_distances, &weighted_distances] (DistanceNodeValue &,
                    DistanceNodeValue &,
                    unsigned long u,
                    double w) {
                unweighted_distances.insert(u);
                weighted_distances.insert(w);
            };
            this->calc_pairwise_tip_distances(f);
            this->unweighted_pairwise_tip_distance_profile_.set_data(unweighted_distances.begin(),
                    unweighted_distances.end());
            this->weighted_pairwise_tip_distance_profile_.set_data(weighted_distances.begin(),
                    weighted_distances.end());
        }

        // This calculates the distances, but does not actually store them.
        // Client code eshould pass in appropriate storage behavior by using a
        // functor for the parameter `store_pairwise_tip_dist_fn`:
        //
        //      f(DistanceNodeValue & t1, DistanceNodeValue & t2, unsigned long steps, double weighted_steps)
        //
        // where first two arguments are the two tips, followed by the number
        // of edges connecting the two, followed by the sum of edge weights
        // along the edges.
        template <typename T>
        void calc_pairwise_tip_distances(T store_pairwise_tip_dist_fn) {
            this->total_tree_length_ = 0.0;
            this->number_of_tips_ = 0;
            // std::multiset<unsigned long> unweighted_distances;
            // std::multiset<double> weighted_distances;
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
                                    store_pairwise_tip_dist_fn(*desc1_nd, *desc2_nd, steps, weight);
                                    // this->weighted_pairwise_tip_distance_[std::make_pair(&(*desc1_nd), &(*desc2_nd))] = weight;
                                    // this->unweighted_pairwise_tip_distance_[std::make_pair(&(*desc1_nd), &(*desc2_nd))] = steps;
                                    // unweighted_distances.insert(steps);
                                    // weighted_distances.insert(weight);
                                }
                            }
                        }
                    }
                }
            }
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

            // this->profile_.add_subprofile_source(this->UNWEIGHTED_PAIRWISE_TIP,
            //         unweighted_distances.begin(),
            //         unweighted_distances.end(),
            //         unweighted_distances.size(),
            //         false);
            // this->profile_.add_subprofile_source(this->WEIGHTED_PAIRWISE_TIP,
            //         weighted_distances.begin(),
            //         weighted_distances.end(),
            //         weighted_distances.size(),
            //         false);

    private:
        Profile                 unweighted_pairwise_tip_distance_profile_;
        Profile                 weighted_pairwise_tip_distance_profile_;
        unsigned long           number_of_tips_;
        double                  total_tree_length_;
        SizesSetType            subtree_leaf_set_sizes_;

}; // DistanceTree

} // namespace pstrudel
#endif
