#ifndef PSTRUDEL_CHARACTER_HPP
#define PSTRUDEL_CHARACTER_HPP

#include <array>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <colugo/console.hpp>
#include <colugo/assert.hpp>

namespace pstrudel {

//////////////////////////////////////////////////////////////////////////////
// Typedefs

typedef int CharacterStateType;
typedef std::vector<CharacterStateType> CharacterStateVectorType;

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequence

class NucleotideSequence {

    public:
        NucleotideSequence();
        NucleotideSequence(const std::string& label);
        ~NucleotideSequence();
        void reserve(unsigned long size) {
            this->sequence_.reserve(size);
        }
        inline unsigned long size() const {
            return this->sequence_.size();
        }
        inline CharacterStateVectorType::iterator begin() {
            return this->sequence_.begin();
        }
        inline CharacterStateVectorType::iterator end() {
            return this->sequence_.end();
        }
        inline const CharacterStateVectorType::const_iterator cbegin() const {
            return this->sequence_.cbegin();
        }
        inline const CharacterStateVectorType::const_iterator cend() const {
            return this->sequence_.cend();
        }
        inline void append_state(CharacterStateType state) {
            this->sequence_.push_back(state);
            auto state_partials = this->state_to_partials_map_.find(state)->second;
            this->partials_.insert(this->partials_.end(), state_partials.begin(), state_partials.end());
        }
        inline void append_state_by_symbol(char s) {
            auto state = NucleotideSequence::get_state_from_symbol(s);
            this->append_state(state);
        }
        inline void append_states_by_symbols(const std::string& s) {
            this->sequence_.reserve(this->sequence_.size() + s.size());
            for (auto & c : s) {
                this->append_state_by_symbol(c);
            }
        }
        inline const CharacterStateType * state_data() const {
            return this->sequence_.data();
        }
        inline const std::string& get_label() const {
            return this->label_;
        }
        inline void set_label(const std::string& label) {
            this->label_ = label;
        }

        void write_states_as_symbols(std::ostream& out) const;
        void write_states_as_symbols(std::ostream& out,
                const CharacterStateVectorType::const_iterator& begin,
                const CharacterStateVectorType::const_iterator& end) const;
        std::string get_states_as_symbols() const;

    protected:
        std::string                 label_;
        CharacterStateVectorType    sequence_;
        std::vector<double>         partials_;

    public:
        static const std::map<char, CharacterStateType>                     symbol_to_state_map_;
        static const std::map<CharacterStateType, char>                     state_to_symbol_map_;
        static const std::map<CharacterStateType, std::array<double, 4>>    state_to_partials_map_;
        static const CharacterStateType                                     missing_data_state;
        static const std::array<double, 4>                                  missing_data_partials;

    public:
        inline static const CharacterStateType get_state_from_symbol(char s) {
            auto state_lookup = NucleotideSequence::symbol_to_state_map_.find(s);
            if (state_lookup == NucleotideSequence::symbol_to_state_map_.end()) {
                colugo::console::abort("Invalid state symbol '", s, "'");
                return 4;
            } else {
                return state_lookup->second;
            }
        }
        inline static const char get_symbol_from_state(CharacterStateType s) {
            auto symbol_lookup = NucleotideSequence::state_to_symbol_map_.find(s);
            if (symbol_lookup == NucleotideSequence::state_to_symbol_map_.end()) {
                colugo::console::abort("Invalid state: ", s);
            }
            return symbol_lookup->second;
        }
        static std::vector<NucleotideSequence> read_fasta(std::istream& src);

}; // NucleotideSequence

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequences

class NucleotideSequences {

    public:
        NucleotideSequences();
        ~NucleotideSequences();
        void clear();
        inline std::vector<NucleotideSequence *>::iterator begin() {
            return this->sequences_.begin();
        }
        inline std::vector<NucleotideSequence *>::iterator end() {
            return this->sequences_.end();
        }
        inline const std::vector<NucleotideSequence *>::const_iterator cbegin() const {
            return this->sequences_.cbegin();
        }
        inline const std::vector<NucleotideSequence *>::const_iterator cend() const {
            return this->sequences_.cend();
        }
        inline unsigned long size() const {
            return this->sequences_.size();
        }
        inline NucleotideSequence * new_sequence(const std::string& label) {
            NucleotideSequence * v = new NucleotideSequence(label);
            this->sequences_.push_back(v);
            this->label_sequence_map_[label] = v;
            return v;
        }
        inline NucleotideSequence * get_sequence(unsigned long index) const {
            COLUGO_ASSERT(index < this->sequences_.size());
            return this->sequences_[index];
        }
        inline NucleotideSequence * get_sequence(const std::string& label) const {
            auto label_sequence = this->label_sequence_map_.find(label);
            COLUGO_ASSERT(label_sequence != this->label_sequence_map_.end());
            return label_sequence->second;
        }
        inline unsigned long get_num_sequences() {
            return this->sequences_.size();
        }
        inline unsigned long get_num_sites() {
            if (this->sequences_.size() > 0) {
                return this->sequences_[0]->size();
            } else {
                return 0;
            }
        }
        void read_fasta(std::istream& src);

    protected:
        std::vector<NucleotideSequence *>               sequences_;
        std::map<std::string, NucleotideSequence *>     label_sequence_map_;

}; // NucleotideSequences

} // namespace pstrudel

#endif
