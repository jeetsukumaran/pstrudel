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
#include <platypus/model/datatable.hpp>
#include <platypus/model/standardinterface.hpp>
#include <platypus/utility/stream.hpp>
#include "profile.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// DistanceNodeValue

class DistanceNodeValue : public platypus::StandardNodeValue<double> {
    public:
        DistanceNodeValue()
            : num_leaves_(0)
            , root_distance_(0.0)
            , age_(0.0) {
        }
        DistanceNodeValue(DistanceNodeValue&& other)
            : platypus::StandardNodeValue<double>(std::move(other))
            , desc_path_lens_(std::move(other.desc_path_lens_))
            , num_leaves_(other.num_leaves_)
            , root_distance_(other.root_distance_)
            , age_(other.age_) {
        }
        DistanceNodeValue(const DistanceNodeValue& other)
            : platypus::StandardNodeValue<double>(other)
            , desc_path_lens_(other.desc_path_lens_)
            , num_leaves_(other.num_leaves_)
            , root_distance_(other.root_distance_)
            , age_(other.age_) {
        }
        DistanceNodeValue& operator=(const DistanceNodeValue& other) {
            platypus::StandardNodeValue<double>::operator=(other);
            this->desc_path_lens_ = other.desc_path_lens_;
            this->num_leaves_ = other.num_leaves_;
            this->root_distance_ = other.root_distance_;
            this->age_ = 0.0;
            return *this;
        }
        void clear() override {
            platypus::StandardNodeValue<double>::clear();
            this->desc_path_lens_.clear();
            this->num_leaves_ = 0;
            this->root_distance_ = 0.0;
            this->age_ = 0.0;
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
        void set_root_distance(double d) {
            this->root_distance_ = d;
        }
        double get_root_distance() const {
            return this->root_distance_;
        }
        void set_age(double d) {
            this->age_ = d;
        }
        double get_age() const {
            return this->age_;
        }
    private:
        std::map<DistanceNodeValue *, std::pair<double, unsigned long>>    desc_path_lens_;
        unsigned long                                                      num_leaves_;
        double                                                             root_distance_;
        double                                                             age_;
}; // DistanceNodeValue

//////////////////////////////////////////////////////////////////////////////
// forward declaration

class DistanceTree;

//////////////////////////////////////////////////////////////////////////////
// PairwiseTipDistanceProfileCalculator

class PairwiseTipDistanceProfileCalculator {

    public:
        PairwiseTipDistanceProfileCalculator(DistanceTree & tree)
            : tree_(tree) { }
        PairwiseTipDistanceProfileCalculator & operator=(const PairwiseTipDistanceProfileCalculator & other);
        double get_unweighted_distance(PairwiseTipDistanceProfileCalculator & other);
        double get_weighted_distance(PairwiseTipDistanceProfileCalculator & other);
        double get_scaled_weighted_distance(PairwiseTipDistanceProfileCalculator & other);

    private:
        void build_pairwise_tip_distance_profiles();

    private:
        DistanceTree &    tree_;
        Profile           unweighted_pairwise_tip_distance_profile_;
        Profile           weighted_pairwise_tip_distance_profile_;
        Profile           scaled_weighted_pairwise_tip_distance_profile_;

}; // PairwiseTipDistanceProfileCalculator

//////////////////////////////////////////////////////////////////////////////
// SymmetricDifferenceCalculator

class SymmetricDifferenceCalculator {
    public:
        typedef std::unordered_multiset<unsigned long> SizesSetType;

    public:
        SymmetricDifferenceCalculator(DistanceTree & tree)
            : tree_(tree) { }
        SymmetricDifferenceCalculator & operator=(const SymmetricDifferenceCalculator & other);
        void calc_subtree_leaf_set_sizes();
        unsigned long calc_leaf_set_sizes_unlabeled_symmetric_difference(SymmetricDifferenceCalculator & other);
        unsigned long get_unlabeled_symmetric_difference(SymmetricDifferenceCalculator & other);
        unsigned long get_unweighted_labeled_symmetric_difference(SymmetricDifferenceCalculator & other);
        unsigned long get_weighted_labeled_symmetric_difference(SymmetricDifferenceCalculator & other);
        static unsigned long calc_set_symmetric_difference(
                const SizesSetType & set1,
                const SizesSetType & set2,
                SizesSetType * common = nullptr,
                SizesSetType * unmatched1 = nullptr,
                SizesSetType * unmatched2 = nullptr);

    private:
        DistanceTree &    tree_;
        SizesSetType      subtree_leaf_set_sizes_;
}; // SymmetricDifferenceCalculator

//////////////////////////////////////////////////////////////////////////////
// LineageThroughTimeProfileCalculator

class LineageThroughTimeProfileCalculator {

    public:
        LineageThroughTimeProfileCalculator(DistanceTree & tree)
            : tree_(tree)
            , lineage_accumulation_through_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , lineage_splitting_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , scaled_lineage_splitting_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , max_leaf_distance_(0.0)
            , min_edge_length_(-1.0) { }
        LineageThroughTimeProfileCalculator & operator=(const LineageThroughTimeProfileCalculator & other);
        double get_lineage_accumulation_profile_distance(LineageThroughTimeProfileCalculator & other);
        double get_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other);
        double get_scaled_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other);
        Profile & get_lineage_accumulation_through_time_profile() {
            return this->lineage_accumulation_through_time_profile_;
        }
        Profile & get_lineage_splitting_time_profile() {
            return this->lineage_splitting_time_profile_;
        }
        Profile & get_scaled_lineage_splitting_time_profile() {
            return this->scaled_lineage_splitting_time_profile_;
        }

    public:
        void clear();

    public:
        // internal, but exposed publically for testing
        void calc_node_root_distances();
        unsigned long get_default_num_transects();
        std::vector<double> build_transect_offsets(unsigned long num_transects=0);
        const Profile & build_lineage_accumulation_through_time_profile(const std::vector<double> & transect_offsets);
        const Profile & build_lineage_accumulation_through_time_profile(unsigned long num_transects=0);
        std::pair<const Profile &, const Profile &> build_lineage_splitting_time_profile();

    protected:
        DistanceTree &    tree_;
        Profile           lineage_accumulation_through_time_profile_;
        Profile           lineage_splitting_time_profile_;
        Profile           scaled_lineage_splitting_time_profile_;
        double            max_leaf_distance_;
        double            min_edge_length_;

}; // LineageThroughTimeProfileCalculator


////////////////////////////////////////////////////////////////////////////////
// DistanceTree

class DistanceTree : public platypus::StandardTree<DistanceNodeValue> {

    public:
        typedef DistanceNodeValue value_type;
        typedef typename StandardTree<DistanceNodeValue>::node_type node_type;

    public:

        /////////////////////////////////////////////////////////////////////////
        // Lifecycle and assignment

        DistanceTree(bool is_rooted=true);
        DistanceTree(DistanceTree && other);
        DistanceTree(const DistanceTree & other);
        DistanceTree & operator=(const DistanceTree & other);
        ~DistanceTree() override {}

        /////////////////////////////////////////////////////////////////////////
        // Basic metrics

        void set_num_tips(unsigned long n) {
            this->number_of_tips_ = n;
        }
        unsigned long get_num_tips() const {
            return this->number_of_tips_;
        }
        unsigned long calc_num_tips() {
            this->number_of_tips_ = 0;
            for (auto ndi = this->leaf_begin(); ndi != this->leaf_end(); ++ndi) {
                ++this->number_of_tips_;
            }
            return this->number_of_tips_;
        }
        void set_total_tree_length(double v) {
            this->total_tree_length_ = v;
        }
        double get_total_tree_length() const {
            return this->total_tree_length_;
        }
        double calc_total_tree_length() {
            this->total_tree_length_ = 0.0;
            for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
                this->total_tree_length_ += ndi->get_edge_length();
            }
        }
        std::vector<double> calc_node_ages(bool include_leaves=true);

        /////////////////////////////////////////////////////////////////////////
        // Unary metrics
        double get_pybus_harvey_gamma();
        double get_N_bar();
        double get_colless_tree_imbalance();
        double get_B1();
        double get_treeness();

        /////////////////////////////////////////////////////////////////////////
        // Manipulators

        // type: 0 = mean coalescent, 1 = converse-coalescent, 2 = uniform
        void add_edge_lengths(int regime=0);

        /////////////////////////////////////////////////////////////////////////
        // Calculators

        // pairwise tup profile
        double get_unweighted_pairwise_tip_profile_distance(DistanceTree & other) {
            return this->pairwise_tip_distance_profile_calculator_.get_unweighted_distance(other.pairwise_tip_distance_profile_calculator_);
        }
        double get_weighted_pairwise_tip_profile_distance(DistanceTree & other) {
            return this->pairwise_tip_distance_profile_calculator_.get_weighted_distance(other.pairwise_tip_distance_profile_calculator_);
        }
        double get_scaled_weighted_pairwise_tip_profile_distance(DistanceTree & other) {
            return this->pairwise_tip_distance_profile_calculator_.get_scaled_weighted_distance(other.pairwise_tip_distance_profile_calculator_);
        }

        // lineage through time
        LineageThroughTimeProfileCalculator & get_lineage_through_time_calculator() {
            return this->lineage_through_time_calculator_;
        }
        double get_lineage_accumulation_profile_distance(DistanceTree & other) {
            return this->lineage_through_time_calculator_.get_lineage_accumulation_profile_distance(other.lineage_through_time_calculator_);
        }
        double get_lineage_splitting_time_profile_distance(DistanceTree & other) {
            return this->lineage_through_time_calculator_.get_lineage_splitting_time_profile_distance(other.lineage_through_time_calculator_);
        }
        double get_scaled_lineage_splitting_time_profile_distance(DistanceTree & other) {
            return this->lineage_through_time_calculator_.get_scaled_lineage_splitting_time_profile_distance(other.lineage_through_time_calculator_);
        }

        // symmetric difference
        SymmetricDifferenceCalculator & get_symmetric_difference_calculator() {
            return this->symmetric_difference_calculator_;
        }
        unsigned long get_unlabeled_symmetric_difference(DistanceTree & other) {
            return this->symmetric_difference_calculator_.get_unlabeled_symmetric_difference(other.symmetric_difference_calculator_);
        }

        template <class R>
        void tabulate_unary_statistics(const std::string & prefix, R & row) {
            row.set(prefix + "ntips"             , this->number_of_tips_);
            row.set(prefix + "length"            , this->total_tree_length_);
            row.set(prefix + "B1"                , this->get_B1());
            row.set(prefix + "colless.imbalance" , this->get_colless_tree_imbalance());
            row.set(prefix + "gamma"             , this->get_pybus_harvey_gamma());
            row.set(prefix + "N.bar"             , this->get_N_bar());
            row.set(prefix + "treeness"          , this->get_treeness());
        }

        template <class T, class R>
        void tabulate_distances(
                T & other_tree,
                R & row,
                bool scale_by_tree_length,
                bool calculate_symmetric_diff) {
            double d = 0.0;
            d = this->get_unweighted_pairwise_tip_profile_distance(other_tree);
            row.set("pwtd.uw", d);
            d = this->get_lineage_accumulation_profile_distance(other_tree);
            row.set("ltt", d);
            if (scale_by_tree_length) {
                d = this->get_scaled_weighted_pairwise_tip_profile_distance(other_tree);
                row.set("pwtd", d);
                d = this->get_scaled_lineage_splitting_time_profile_distance(other_tree);
                row.set("lst", d);
            } else {
                d = this->get_weighted_pairwise_tip_profile_distance(other_tree);
                row.set("pwtd", d);
                d = this->get_lineage_splitting_time_profile_distance(other_tree);
                row.set("lst", d);
            }
            if (calculate_symmetric_diff) {
                d = this->get_unlabeled_symmetric_difference(other_tree);
                row.set("rfdu", d);
            }
        }

    public:

        static void add_unary_statistic_columns(
                const std::string & prefix,
                platypus::DataTable & table,
                platypus::stream::OutputStreamFormatters & col_formatters,
                bool is_key_column) {
            table.add_column<unsigned long>(prefix + "ntips", {}, is_key_column);
            table.add_column<double>(prefix + "length", col_formatters, is_key_column);
            table.add_column<double>(prefix + "B1", col_formatters, is_key_column);
            table.add_column<double>(prefix + "colless.imbalance", col_formatters, is_key_column);
            table.add_column<double>(prefix + "gamma", col_formatters, is_key_column);
            table.add_column<double>(prefix + "N.bar", col_formatters, is_key_column);
            table.add_column<double>(prefix + "treeness", col_formatters, is_key_column);
        }

        static void add_results_data_columns(
                platypus::DataTable & table,
                platypus::stream::OutputStreamFormatters & col_formatters,
                bool calculate_symmetric_diff) {
            for (auto & y_distance_name : DistanceTree::tree_pattern_y_distance_names_) {
                table.add_data_column<double>(y_distance_name, col_formatters);
            }
            if (calculate_symmetric_diff) {
                table.add_data_column<double>("rfdu", col_formatters);
            }
        }

        /////////////////////////////////////////////////////////////////////////
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
            this->number_of_tips_ = 0;
            this->total_tree_length_ = 0.0;
            for (auto ndi = this->postorder_begin(); ndi != this->postorder_end(); ++ndi) {
                this->total_tree_length_ += ndi->get_edge_length();
                if (ndi.is_leaf()) {
                    ndi->set_desc_path_len(*ndi, 0.0, 0);
                    ++this->number_of_tips_;
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
                                }
                            }
                        }
                    }
                }
            }
        }

    private:
        std::vector<DistanceTree::node_type *> get_nodes_in_level_order();

    private:
        unsigned long                            number_of_tips_;
        double                                   total_tree_length_;
        PairwiseTipDistanceProfileCalculator     pairwise_tip_distance_profile_calculator_;
        SymmetricDifferenceCalculator            symmetric_difference_calculator_;
        LineageThroughTimeProfileCalculator      lineage_through_time_calculator_;
        static const std::vector<std::string>    tree_pattern_y_distance_names_;
        double                                   B1_;
        double                                   colless_tree_imbalance_;
        double                                   pybus_harvey_gamma_;
        double                                   N_bar_;
        double                                   treeness_;

}; // DistanceTree



} // namespace pstrudel
#endif
