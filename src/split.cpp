#include <algorithm>
#include "split.hpp"

namespace pstrudel {

////////////////////////////////////////////////////////////////////////////////
// BitMask

BitMask::BitMask() {
}

BitMask::BitMask(const BitMask & other) {
    this->ulong_val_ = other.ulong_val_;
    this->bits_ = other.bits_;
}

BitMask::BitMask(unsigned long ulong_val) {
    this->set_from_ulong(ulong_val);
}

void BitMask::set(unsigned long ulong_val) {
    this->set_from_ulong(ulong_val);
}

void BitMask::set_from_ulong(unsigned long ulong_val) {
    this->bits_.clear();
    while (ulong_val) {
        if (ulong_val & 1)
            this->bits_.push_back(true);
        else
            this->bits_.push_back(false);
        ulong_val >>= 1;
    }
    std::reverse(this->bits_.begin(), this->bits_.end());
}

void BitMask::set_for_index(unsigned long idx) {
    this->bits_.clear();
    this->bits_.push_back(true);
    for (unsigned long i = 0; i < idx-1; ++i) {
        this->bits_.push_back(false);
    }
}

////////////////////////////////////////////////////////////////////////////////
// Split

Split::Split() {
}

Split::Split(const BitMask& split_hash, double edge_length)
    : split_hash_(split_hash)
      , edge_length_(edge_length) {
      }

Split::Split(const Split& other)
    : split_hash_(other.split_hash_)
      , edge_length_(other.edge_length_) {
}

Split::Split(Split&& other)
    : split_hash_(std::move(other.split_hash_))
      , edge_length_(std::move(other.edge_length_)) {
}

Split::~Split() {
}

Split & Split::operator=(const Split& other) {
    this->split_hash_ = other.split_hash_;
    this->edge_length_ = other.edge_length_;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
// TaxonSet

SplitHashType & TaxonSet::add_taxon(const std::string& taxon_label) {
    std::string normalized = taxon_label;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), ::toupper);
    auto found = this->taxon_map_.find(normalized);
    if (found == this->taxon_map_.end()) {
        unsigned long idx = this->taxon_map_.size();
        auto * hash_ptr = new SplitHashType;
        hash_ptr->set_for_index(idx+1);
        this->taxon_map_[normalized] = hash_ptr;
        return *hash_ptr;
    } else {
        return *(found->second);
    }
}

SplitHashType & TaxonSet::get_taxon_hash(const std::string& taxon_label) {
    return this->add_taxon(taxon_label);
}

} // namespace pstrudel

