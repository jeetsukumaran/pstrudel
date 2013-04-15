#ifndef PSTRUDEL_SYMMETRIC_DISTANCE_TREE_HPP
#define PSTRUDEL_SYMMETRIC_DISTANCE_TREE_HPP

#include <unordered_set>
#include "basic_tree.hpp"

namespace pstrudel {

///////////////////////////////////////////////////////////////////////////////
// SymmetricDifferenceNodeValue

class SymmetricDifferenceNodeValue : public BasicNodeValue {
    public:
        SymmetricDifferenceNodeValue()
            : num_leaves_(0)
              , num_descendent_nodes_(0) {
        }
        SymmetricDifferenceNodeValue(const SymmetricDifferenceNodeValue & other)
            : num_leaves_(other.num_leaves_)
              , num_descendent_nodes_(other.num_descendent_nodes_) {
        }
        template <class T> SymmetricDifferenceNodeValue & operator=(const T & other) {
            return *this;
        }
        unsigned long get_num_leaves() const {
            return this->num_leaves_;
        }
        void set_num_leaves(unsigned long num_leaves) {
            this->num_leaves_ = num_leaves;
        }
        unsigned long get_num_descendent_nodes() const {
            return this->num_descendent_nodes_;
        }
        void set_num_descendent_nodes(unsigned long num_descendent_nodes) {
            this->num_descendent_nodes_ = num_descendent_nodes;
        }
    private:
        unsigned long num_leaves_;
        unsigned long num_descendent_nodes_;
}; // SymmetricDifferenceNodeValue

///////////////////////////////////////////////////////////////////////////////
// SymmetricDifferenceTree

class SymmetricDifferenceTree : public BasicTree<SymmetricDifferenceNodeValue> {

    public:
        typedef std::unordered_multiset<unsigned long> SizesSetType;

    public:
        SymmetricDifferenceTree();
        SymmetricDifferenceTree(const SymmetricDifferenceTree & other);
        template <class T> SymmetricDifferenceTree(const T & other)
            : BasicTree<SymmetricDifferenceNodeValue>(other) {
        }
        void calc_subtree_sizes();
        unsigned long calc_leaf_set_sizes_unlabeled_symmetric_difference(SymmetricDifferenceTree & other);
        unsigned long calc_clade_sizes_unlabeled_symmetric_difference(SymmetricDifferenceTree & other);
        unsigned long get_unlabeled_symmetric_difference(SymmetricDifferenceTree & other);
        unsigned long get_unweighted_labeled_symmetric_difference(SymmetricDifferenceTree & other);
        unsigned long get_weighted_labeled_symmetric_difference(SymmetricDifferenceTree & other);

    public:
        static unsigned long calc_set_symmetric_difference(
                const SizesSetType & set1,
                const SizesSetType & set2,
                SizesSetType * common=nullptr,
                SizesSetType * unmatched1=nullptr,
                SizesSetType * unmatched2=nullptr);

    private:
        SizesSetType        subtree_leaf_set_sizes_;
        SizesSetType        subtree_clade_sizes_;
}; // SymmetricDifferenceTree

} // namespace pstrudel

#endif
