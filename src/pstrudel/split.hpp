#ifndef PSTRUDEL_SPLIT_HPP
#define PSTRUDEL_SPLIT_HPP

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <bitset>
#include <colugo/utility.hpp>

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// BitMask

class BitMask {

    public:
        typedef std::vector<bool> BitVectorType;

    public:
        BitMask();
        BitMask(const BitMask & other);
        BitMask(unsigned long ulong_val);
        void set(unsigned long ulong_val);
        void set_from_ulong(unsigned long ulong_val);
        void set_for_index(unsigned long idx);
        BitMask& operator=(const BitMask & other) {
            this->ulong_val_ = other.ulong_val_;
            this->bits_ = other.bits_;
            return *this;
        }
        BitMask& operator|=(const BitMask& rhs) {
            unsigned long s1 = 0;
            unsigned long s2 = 0;
            if (this->bits_.size() > rhs.bits_.size()) {
                s1 = this->bits_.size() - rhs.bits_.size();
            } else if (rhs.bits_.size() > this->bits_.size()) {
                unsigned long x = rhs.bits_.size() - this->bits_.size();
                this->bits_.insert(this->bits_.begin(), x, false);
                COLUGO_ASSERT(this->bits_.size() == rhs.bits_.size());
                // this->bits_.insert(this->bits_.begin(), rhs.bits_.begin(), rhs.bits_.begin() + x);
                // s1 = x;
                // s2 = x;
            }
            while (s1 < this->bits_.size() && s2 < rhs.bits_.size()) {
                this->bits_[s1] = this->bits_[s1] || rhs.bits_[s2];
                ++s1;
                ++s2;
            }
            // std::transform(this->bits_.begin(), this->bits_.end(),
            //    rhs.bits_.begin(), this->bits_.begin(), std::logical_or<bool>());
            return *this;
        }
        bool operator<(const BitMask& rhs) const {
            // return this->ulong_val_ < rhs.ulong_val_;
            return this->bits_ < rhs.bits_;
        }
        bool operator==(const BitMask& rhs) const {
            // return this->ulong_val_ == rhs.ulong_val_;
            return this->bits_ == rhs.bits_;
        }
        void dump(std::ostream& out) const {
            for (BitVectorType::const_iterator bi = this->bits_.begin();
                    bi != this->bits_.end();
                    ++ bi) {
                out << *bi;
            }
        }

    private:
        unsigned long       ulong_val_;
        BitVectorType       bits_;
        std::hash<BitVectorType>      bitmask_hash_f_;

}; // BitMask

////////////////////////////////////////////////////////////////////////////////
// Split

class Split {

    public:
        typedef BitMask SplitHashType;

    public:
        Split();
        Split(const BitMask& split_hash, double edge_length=0.0);
        Split(const Split& other);
        Split(Split&& other);
        ~Split();
        Split& operator=(const Split& other);
        inline std::size_t bitmask_hash() {
            return 0;
        }

    private:
        BitMask                 split_hash_;
        double                  edge_length_;

}; // Split

typedef BitMask SplitHashType;
typedef std::map<SplitHashType, double> TreeSplitMapType;

////////////////////////////////////////////////////////////////////////////////
// TaxonSet

class TaxonSet {

    public:
        typedef Split::SplitHashType SplitHashType;

    public:
        TaxonSet() { }
        ~TaxonSet() {
            for (auto & tm : this->taxon_map_) {
                if (tm.second != nullptr) {
                    delete tm.second;
                }
            }
        }
        SplitHashType & add_taxon(const std::string& taxon_label);
        SplitHashType & get_taxon_hash(const std::string& taxon_label);

    private:
        std::map<std::string, SplitHashType *>   taxon_map_;

}; // TaxonSet

////////////////////////////////////////////////////////////////////////////////
// TreeComparisonResult
struct TreeComparisonResult {
    public:
        unsigned long   false_positives;
        unsigned long   false_negatives;
        unsigned long   symmetric_difference;
        double          edge_len_euclidean_distance;
        double          weighted_rf_distance;
}; // TreeComparisonResult

////////////////////////////////////////////////////////////////////////////////
// TreeComparisonCalculator

template <class T>
class TreeComparisonCalculator {


    public:
        typedef T                               TreeType;
        typedef typename TreeType::node_type    NodeType;
        typedef typename TreeType::value_type   NodeValueType;
        typedef Split::SplitHashType            SplitHashType;
        typedef std::function<void (NodeValueType &, const SplitHashType &)> SetNodeSplitFuncType;

    public:

        TreeComparisonCalculator() {
        }

        TreeComparisonCalculator(SetNodeSplitFuncType set_node_split_func)
            : set_node_split_func_(set_node_split_func) {
        }

        TreeSplitMapType& calc_tree_splits(const TreeType & tree) {
            auto & current_split_map = this->tree_splits_[&tree];
            std::map<typename TreeType::value_type *, SplitHashType> node_split_map;
            for (auto ndi = tree.postorder_begin(); ndi != tree.postorder_end(); ++ndi) {
                if (ndi.is_leaf()) {
                    COLUGO_ASSERT(node_split_map.find(&(*ndi)) == node_split_map.end());
                    node_split_map[&(*ndi)] = this->taxa_.get_taxon_hash(ndi->get_label());
                    auto & hash = node_split_map[&(*ndi)];
                    COLUGO_ASSERT(current_split_map.find(hash) == current_split_map.end());
                    current_split_map.insert(std::make_pair(hash, ndi->get_edge_length()));
                    if (this->set_node_split_func_) {
                        this->set_node_split_func_(*ndi, hash);
                    }
                } else {
                    auto chi = tree.children_begin(ndi);
                    COLUGO_ASSERT(node_split_map.find(&(*chi)) != node_split_map.end());
                    COLUGO_ASSERT(node_split_map.find(&(*ndi)) == node_split_map.end());
                    node_split_map[&(*ndi)] = node_split_map[&(*chi)];      // initialize with bitmask of first child
                    auto & hash = node_split_map[&(*ndi)];
                    for (++chi; chi != tree.children_end(ndi); ++chi) {
                        hash |= node_split_map[&(*chi)];                    // OR together remaining children
                    }
                    COLUGO_ASSERT(current_split_map.find(hash) == current_split_map.end());
                    current_split_map.insert(std::make_pair(hash, ndi->get_edge_length()));
                    if (this->set_node_split_func_) {
                        this->set_node_split_func_(*ndi, hash);
                    }
                }
            }
            return current_split_map;
        }

        TreeSplitMapType& get_tree_split_map(const TreeType& tree) {
            auto tree_split_map_iter = this->tree_splits_.find(&tree);
            if (tree_split_map_iter == this->tree_splits_.end()) {
                return this->calc_tree_splits(tree);
            }
            return this->tree_splits_[&tree];
        }

        TreeComparisonResult compare_trees(const TreeType& tree1, const TreeType& tree2) {
            TreeComparisonResult result;
            auto & split_set1 = this->get_tree_split_map(tree1);
            auto & split_set2 = this->get_tree_split_map(tree2);
            result.false_positives = 0;
            result.false_negatives = 0;
            result.symmetric_difference = 0;
            result.edge_len_euclidean_distance = 0.0;
            result.weighted_rf_distance = 0.0;
            for (auto & sp1_iter : split_set1) {
                auto & split1 = sp1_iter.first;
                auto & edge1 = sp1_iter.second;
                auto sp2_iter = split_set2.find(split1);
                if (sp2_iter == split_set2.end()) {
                    result.false_negatives += 1;
                    result.edge_len_euclidean_distance += std::pow(edge1, 2);
                    result.weighted_rf_distance += edge1;
                } else {
                    auto & edge2 = sp2_iter->second;
                    double v = edge1 - edge2;
                    result.edge_len_euclidean_distance += std::pow(v, 2);
                    result.weighted_rf_distance += std::abs(v);
                }
            }
            for (auto & sp2_iter : split_set2) {
                auto & split2 = sp2_iter.first;
                auto & edge2 = sp2_iter.second;
                auto sp1_iter = split_set1.find(split2);
                if (sp1_iter == split_set1.end()) {
                    result.false_positives += 1;
                    result.edge_len_euclidean_distance += std::pow(edge2, 2);
                    result.weighted_rf_distance += edge2;
                }
            }
            result.symmetric_difference = result.false_positives + result.false_negatives;
            result.edge_len_euclidean_distance = std::sqrt(result.edge_len_euclidean_distance);
            return result;
        }

    private:
        TaxonSet                                        taxa_;
        std::map<const TreeType *, TreeSplitMapType>    tree_splits_;
        SetNodeSplitFuncType                            set_node_split_func_;

}; // TreeComparisonCalculator


} // namespace pstrudel

#endif
