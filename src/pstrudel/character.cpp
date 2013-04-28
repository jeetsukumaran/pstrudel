#include <algorithm>
#include <iterator>
#include <colugo/utility.hpp>
#include "character.hpp"

namespace pstrudel {

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequence

const std::map<char, CharacterStateType> NucleotideSequence::symbol_to_state_map_ {
    {'A', 0},
    {'C', 1},
    {'G', 2},
    {'T', 3},
    {'U', 3},
    {'N', 4},
    {'X', 4},
    {'-', 4},
    {'?', 4},
    {'R', 5},
    {'Y', 6},
    {'M', 7},
    {'W', 8},
    {'S', 9},
    {'K', 10},
    {'V', 11},
    {'H', 12},
    {'D', 13},
    {'B', 14}
};

const std::map<CharacterStateType, char> NucleotideSequence::state_to_symbol_map_ {
    { 0, 'A'},
    { 1, 'C'},
    { 2, 'G'},
    { 3, 'T'},
    // { 4, 'N'},
    // { 4, 'X'},
    { 4, '-'},
    // { 4, '?'},
    { 5, 'R'},
    { 6, 'Y'},
    { 7, 'M'},
    { 8, 'W'},
    { 9, 'S'},
    {10, 'K'},
    {12, 'H'},
    {13, 'D'},
    {14, 'B'},
    {11, 'V'}
};

const std::map<CharacterStateType, std::array<double, 4>> NucleotideSequence::state_to_partials_map_ {
    { 0, {{1.0, 0.0, 0.0, 0.0}}},   // A
    { 1, {{0.0, 1.0, 0.0, 0.0}}},   // C
    { 2, {{0.0, 0.0, 1.0, 0.0}}},   // G
    { 3, {{0.0, 0.0, 0.0, 1.0}}},   // T, U
    { 4, {{1.0, 1.0, 1.0, 1.0}}},   // N, X, -, ?
    { 5, {{1.0, 0.0, 1.0, 0.0}}},   // R
    { 6, {{0.0, 1.0, 0.0, 1.0}}},   // Y
    { 7, {{1.0, 1.0, 0.0, 0.0}}},   // M
    { 8, {{1.0, 0.0, 0.0, 1.0}}},   // W
    { 9, {{0.0, 1.0, 1.0, 0.0}}},   // S
    {10, {{0.0, 0.0, 1.0, 1.0}}},   // K
    {11, {{1.0, 1.0, 1.0, 0.0}}},   // V
    {12, {{1.0, 1.0, 0.0, 1.0}}},   // H
    {13, {{1.0, 0.0, 1.0, 1.0}}},   // D
    {14, {{0.0, 1.0, 1.0, 1.0}}}    // B
};

const CharacterStateType NucleotideSequence::missing_data_state(NucleotideSequence::symbol_to_state_map_.find('?')->second);
const std::array<double, 4> NucleotideSequence::missing_data_partials = NucleotideSequence::state_to_partials_map_.find(NucleotideSequence::missing_data_state)->second;

NucleotideSequence::NucleotideSequence() {
}

NucleotideSequence::NucleotideSequence(const std::string& label) :
    label_(label) {
}

NucleotideSequence::~NucleotideSequence() {
}

void NucleotideSequence::write_states_as_symbols(std::ostream& out) const {
    for (auto & s : this->sequence_) {
        out << this->state_to_symbol_map_.find(s)->second;
    }
}

void NucleotideSequence::write_states_as_symbols(std::ostream& out,
    const CharacterStateVectorType::const_iterator& begin,
    const CharacterStateVectorType::const_iterator& end) const {
    for (auto iter = begin; iter != end; ++iter) {
        out << this->state_to_symbol_map_.find(*iter)->second;
    }
}

std::string NucleotideSequence::get_states_as_symbols() const {
    std::ostringstream out;
    this->write_states_as_symbols(out);
    return out.str();
}

// static std::vector<NucleotideSequence> read_fasta(std::istream& src) {
//     unsigned long line_idx = 0;
//     std::vector<NucleotideSequence> seqs;
//     NucleotideSequence * seq = nullptr;
//     for(std::string line; std::getline(src, line); ++line_idx) {
//         if (line.empty()) {
//             continue;
//         }
//         if (line[0] == '>') {
//             seqs.emplace_back(line.substr(1, line.size()));
//             seq = &(seqs.back());
//         } else {
//             for (auto & c : line) {
//                 if (seq == nullptr) {
//                     colugo::colugo_abort("FASTA file read error: Line ",
//                             line_idx+1,
//                             "Expecting sequence label (i.e., line starting with '>')");
//                 }
//                 if (!std::isspace(c)) {
//                     seq->append_state_by_symbol(c);
//                 }
//             }
//         }
//     }
//     return seqs;
// }

//////////////////////////////////////////////////////////////////////////////
// NucleotideSequences

NucleotideSequences::NucleotideSequences() {
}

NucleotideSequences::~NucleotideSequences() {
    this->clear();
}

void NucleotideSequences::clear() {
    for (auto & v : this->sequences_) {
        if (v) {
            delete v;
        }
    }
    this->sequences_.clear();
    this->label_sequence_map_.clear();
}

void NucleotideSequences::read_fasta(std::istream& src) {
    unsigned long line_idx = 0;
    NucleotideSequence * seq = nullptr;
    for(std::string line; std::getline(src, line); ++line_idx) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            seq = this->new_sequence(line.substr(1, line.size()));
        } else {
            for (auto & c : line) {
                if (!seq) {
                    colugo::colugo_abort("Expecting sequence label (i.e., line starting with '>')");
                }
                if (!std::isspace(c)) {
                    seq->append_state_by_symbol(c);
                }
            }
        }
    }
}

} // namespace pstrudel


