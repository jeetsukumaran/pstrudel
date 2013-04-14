#ifndef PSTRUDEL_BASIC_TREE_HPP
#define PSTRUDEL_BASIC_TREE_HPP

#include <string>
#include <vector>
#include "colugo-utilities/src/utility.hpp"
#include "platypus-phyloinformary/src/tree.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// BasicNodeValue

class BasicNodeValue {

    public:
        BasicNodeValue(double edge_length=0.0)
            : edge_length_(edge_length) {
        }
        BasicNodeValue(BasicNodeValue && other)
            : label_(std::move(other.label_))
              , edge_length_(std::move(other.edge_length_)) {
        }
        BasicNodeValue(const BasicNodeValue & other)
            : label_(other.label_)
              , edge_length_(other.edge_length_) {
        }
        BasicNodeValue & operator=(const BasicNodeValue & other) {
            this->label_ = other.label_;
            this->edge_length_ = other.edge_length_;
            return *this;
        }
        virtual ~BasicNodeValue() {
        }
        virtual void clear() {
            this->label_ = "";
            this->edge_length_ = 0.0;
        }
        inline const std::string& get_label() const {
            return this->label_;
        }
        inline void set_label(const std::string& label) {
            this->label_ = label;
        }
        inline double get_edge_length() const {
            return this->edge_length_;
        }
        inline void set_edge_length(double edge_length) {
            this->edge_length_ = edge_length;
        }
    private:
        std::string     label_;
        double          edge_length_;

}; // BasicNodeValue

////////////////////////////////////////////////////////////////////////////////
// BasicTree

template <class T>
class BasicTree : public platypus::Tree<T> {
    public:
        BasicTree(bool is_rooted=true)
            : is_rooted_(is_rooted) {
        }
        BasicTree(BasicTree&& other)
            : platypus::Tree<T>(std::move(other))
              , is_rooted_(std::move(other.is_rooted_))
              , label_(std::move(other.label_)) {
        }
        BasicTree(const BasicTree & other)
            : platypus::Tree<T>(other)
              , is_rooted_(other.is_rooted_)
              , label_(other.label_) {
        }
        BasicTree & operator=(const BasicTree & other) {
            platypus::Tree<T>::operator=(other);
            this->is_rooted_ = other.is_rooted_;
            this->label_ = other.label_;
            return *this;
        }
        virtual ~BasicTree() {
        }
        bool is_rooted() const {
            return this->is_rooted_;
        }
        void set_rooted(bool is_rooted) {
            this->is_rooted_ = is_rooted;
        }
        const std::string& get_label() const {
            return this->label_;
        }
        void set_label(const std::string& label) {
            this->label_ = label;
        }
    private:
        bool            is_rooted_;
        std::string     label_;
}; // BasicTree

} // namespace pstrudel

#endif
