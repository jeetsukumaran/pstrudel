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

class TreeShape;

//////////////////////////////////////////////////////////////////////////////
// PairwiseTipDistanceProfileCalculator

class PairwiseTipDistanceProfileCalculator {

    public:
        PairwiseTipDistanceProfileCalculator(TreeShape & tree)
            : tree_(tree) { }
        PairwiseTipDistanceProfileCalculator & operator=(const PairwiseTipDistanceProfileCalculator & other);
        double get_unweighted_distance(PairwiseTipDistanceProfileCalculator & other, bool weight_values_by_profile_size=false);
        double get_weighted_distance(PairwiseTipDistanceProfileCalculator & other, bool weight_values_by_profile_size=false);
        double get_scaled_weighted_distance(PairwiseTipDistanceProfileCalculator & other, bool weight_values_by_profile_size=false);

    private:
        void build_pairwise_tip_distance_profiles();

    private:
        TreeShape &    tree_;
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
        SymmetricDifferenceCalculator(TreeShape & tree)
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
        TreeShape &    tree_;
        SizesSetType      subtree_leaf_set_sizes_;
}; // SymmetricDifferenceCalculator

//////////////////////////////////////////////////////////////////////////////
// LineageThroughTimeProfileCalculator

class LineageThroughTimeProfileCalculator {

    public:
        LineageThroughTimeProfileCalculator(TreeShape & tree)
            : tree_(tree)
            , lineage_accumulation_through_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , lineage_splitting_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , scaled_lineage_splitting_time_profile_(0, Profile::InterpolationMethod::STAIRCASE)
            , max_leaf_distance_(0.0)
            , min_edge_length_(-1.0) { }
        LineageThroughTimeProfileCalculator & operator=(const LineageThroughTimeProfileCalculator & other);
        double get_lineage_accumulation_profile_distance(LineageThroughTimeProfileCalculator & other, bool weight_values_by_profile_size=false);
        double get_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other, bool weight_values_by_profile_size=false);
        double get_scaled_lineage_splitting_time_profile_distance(LineageThroughTimeProfileCalculator & other, bool weight_values_by_profile_size=false);
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
        TreeShape &    tree_;
        Profile           lineage_accumulation_through_time_profile_;
        Profile           lineage_splitting_time_profile_;
        Profile           scaled_lineage_splitting_time_profile_;
        double            max_leaf_distance_;
        double            min_edge_length_;

}; // LineageThroughTimeProfileCalculator


////////////////////////////////////////////////////////////////////////////////
// TreeShape

class TreeShape : public platypus::StandardTree<DistanceNodeValue> {

    public:
        typedef DistanceNodeValue value_type;
        typedef typename StandardTree<DistanceNodeValue>::node_type node_type;

    public:

        /////////////////////////////////////////////////////////////////////////
        // Lifecycle and assignment

        TreeShape(bool is_rooted=true);
        TreeShape(TreeShape && other);
        TreeShape(const TreeShape & other);
        TreeShape & operator=(const TreeShape & other);
        ~TreeShape() override {}

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
        void create_coalescent_intervals(int regime=0);
        void add_edge_lengths(unsigned long length, bool equalize_root_tip_distances);

        /////////////////////////////////////////////////////////////////////////
        // Calculators

        // pairwise tup profile
        double get_unweighted_pairwise_tip_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->pairwise_tip_distance_profile_calculator_.get_unweighted_distance(other.pairwise_tip_distance_profile_calculator_, weight_values_by_profile_size);
        }
        double get_weighted_pairwise_tip_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->pairwise_tip_distance_profile_calculator_.get_weighted_distance(other.pairwise_tip_distance_profile_calculator_, weight_values_by_profile_size);
        }
        double get_scaled_weighted_pairwise_tip_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->pairwise_tip_distance_profile_calculator_.get_scaled_weighted_distance(other.pairwise_tip_distance_profile_calculator_, weight_values_by_profile_size);
        }

        // lineage through time
        LineageThroughTimeProfileCalculator & get_lineage_through_time_calculator() {
            return this->lineage_through_time_calculator_;
        }
        double get_lineage_accumulation_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->lineage_through_time_calculator_.get_lineage_accumulation_profile_distance(other.lineage_through_time_calculator_, weight_values_by_profile_size);
        }
        double get_lineage_splitting_time_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->lineage_through_time_calculator_.get_lineage_splitting_time_profile_distance(other.lineage_through_time_calculator_, weight_values_by_profile_size);
        }
        double get_scaled_lineage_splitting_time_profile_distance(TreeShape & other, bool weight_values_by_profile_size=false) {
            return this->lineage_through_time_calculator_.get_scaled_lineage_splitting_time_profile_distance(other.lineage_through_time_calculator_, weight_values_by_profile_size);
        }

        // symmetric difference
        SymmetricDifferenceCalculator & get_symmetric_difference_calculator() {
            return this->symmetric_difference_calculator_;
        }
        unsigned long get_unlabeled_symmetric_difference(TreeShape & other) {
            return this->symmetric_difference_calculator_.get_unlabeled_symmetric_difference(other.symmetric_difference_calculator_);
        }

        // TODO: refactor!
        // coalescent intervals
        double get_unscaled_coalescent_interval_profile_distance(TreeShape & other, bool weight_values_by_profile_size) {
            if (this->unscaled_coalescent_interval_profile_.empty()) {
                this->build_coalescent_interval_profile();
            }
            if (other.unscaled_coalescent_interval_profile_.empty()) {
                other.build_coalescent_interval_profile();
            }
            return this->unscaled_coalescent_interval_profile_.get_distance(other.unscaled_coalescent_interval_profile_, weight_values_by_profile_size);
        }

        double get_scaled_coalescent_interval_profile_distance(TreeShape & other, bool weight_values_by_profile_size) {
            if (this->scaled_coalescent_interval_profile_.empty()) {
                this->build_coalescent_interval_profile();
            }
            if (other.scaled_coalescent_interval_profile_.empty()) {
                other.build_coalescent_interval_profile();
            }
            return this->scaled_coalescent_interval_profile_.get_distance(other.scaled_coalescent_interval_profile_, weight_values_by_profile_size);
        }

        void build_coalescent_interval_profile() {
            if (this->coalescent_intervals_.empty()) {
                this->calc_coalescent_intervals();
            }
            this->unscaled_coalescent_interval_profile_.set_data(this->coalescent_intervals_.begin(), this->coalescent_intervals_.end(), true);
            double tree_length = this->get_total_tree_length();
            std::vector<double> scaled_ci;
            scaled_ci.reserve(this->coalescent_intervals_.size());
            for (auto & ci : this->coalescent_intervals_) {
                scaled_ci.push_back(ci / tree_length);
            }
            this->scaled_coalescent_interval_profile_.set_data(scaled_ci.begin(), scaled_ci.end(), true);
        }

        void calc_coalescent_intervals() {
            this->coalescent_intervals_.clear();
            if (this->number_of_tips_ == 0) {
                this->calc_num_tips();
            }
            TreeShape::node_type * node_ptr = nullptr;
            auto node_ages = this->calc_node_ages(false);
            if (node_ages.empty()) {
                throw std::runtime_error("Cannot calculate statistic on empty tree");
            }
            this->coalescent_intervals_.reserve(node_ages.size());
            auto cur_age_iter = node_ages.begin();
            double older = node_ages[0];
            for (unsigned long age_idx = 1; age_idx < node_ages.size(); ++age_idx) {
                this->coalescent_intervals_.push_back(older - node_ages[age_idx]);
                older = node_ages[age_idx];
            }
            this->coalescent_intervals_.push_back(older);
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
                const std::string & prefix,
                T & other_tree,
                R & row,
                bool scale_by_tree_length,
                bool calculate_symmetric_diff,
                bool calculate_unary_statistics_differences,
                bool weight_values_by_profile_size) {
            double d = 0.0;
            if (calculate_unary_statistics_differences) {
                row.set("diff." + prefix + "ntips"             , std::abs(this->number_of_tips_ - other_tree.number_of_tips_));
                row.set("diff." + prefix + "length"            , std::abs(this->total_tree_length_ - other_tree.total_tree_length_));
                row.set("diff." + prefix + "B1"                , std::abs(this->get_B1() - other_tree.get_B1()));
                row.set("diff." + prefix + "colless.imbalance" , std::abs(this->get_colless_tree_imbalance() - other_tree.get_colless_tree_imbalance()));
                row.set("diff." + prefix + "gamma"             , std::abs(this->get_pybus_harvey_gamma() - other_tree.get_pybus_harvey_gamma()));
                row.set("diff." + prefix + "N.bar"             , std::abs(this->get_N_bar() - other_tree.get_N_bar()));
                row.set("diff." + prefix + "treeness"          , std::abs(this->get_treeness() - other_tree.get_treeness()));
            }
            d = this->get_unweighted_pairwise_tip_profile_distance(other_tree, weight_values_by_profile_size);
            row.set(prefix + "pwtd.uw", d);
            d = this->get_lineage_accumulation_profile_distance(other_tree, weight_values_by_profile_size);
            row.set(prefix + "ltt", d);
            if (scale_by_tree_length) {
                d = this->get_scaled_weighted_pairwise_tip_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "pwtd", d);
                d = this->get_scaled_lineage_splitting_time_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "lst", d);
                d = this->get_scaled_coalescent_interval_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "coal.intv", d);
            } else {
                d = this->get_weighted_pairwise_tip_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "pwtd", d);
                d = this->get_lineage_splitting_time_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "lst", d);
                d = this->get_unscaled_coalescent_interval_profile_distance(other_tree, weight_values_by_profile_size);
                row.set(prefix + "coal.intv", d);
            }
            if (calculate_symmetric_diff) {
                d = this->get_unlabeled_symmetric_difference(other_tree);
                row.set(prefix + "rfdu", d);
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
                const std::string & prefix,
                platypus::DataTable & table,
                platypus::stream::OutputStreamFormatters & col_formatters,
                bool calculate_symmetric_diff,
                bool calculate_unary_statistics_differences) {
            if (calculate_unary_statistics_differences) {
                TreeShape::add_unary_statistic_columns("diff." + prefix,
                        table,
                        col_formatters,
                        false);
            }
            for (auto & y_distance_name : TreeShape::tree_pattern_y_distance_names_) {
                table.add_data_column<double>(prefix + y_distance_name, col_formatters);
            }
            if (calculate_symmetric_diff) {
                table.add_data_column<double>(prefix + "rfdu", col_formatters);
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
        unsigned long                            number_of_tips_;
        double                                   total_tree_length_;
        PairwiseTipDistanceProfileCalculator     pairwise_tip_distance_profile_calculator_;
        SymmetricDifferenceCalculator            symmetric_difference_calculator_;
        LineageThroughTimeProfileCalculator      lineage_through_time_calculator_;
        std::vector<double>                      coalescent_intervals_;
        Profile                                  unscaled_coalescent_interval_profile_;
        Profile                                  scaled_coalescent_interval_profile_;
        static const std::vector<std::string>    tree_pattern_y_distance_names_;
        double                                   B1_;
        double                                   colless_tree_imbalance_;
        double                                   pybus_harvey_gamma_;
        double                                   N_bar_;
        double                                   treeness_;

}; // TreeShape



} // namespace pstrudel
#endif
