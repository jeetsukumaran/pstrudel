#ifndef PSTRUDEL_HPP
#define PSTRUDEL_HPP
#include "version.h"

#include <string>

namespace pstrudel {

std::string get_program_source_identification() {
#if defined(PROJECT_SOURCE_IDENTIFIED)
    return PROJECT_SOURCE_IDENTIFIER;
#else
    return std::string(__DATE__) + __TIME__;
#endif
}

std::string get_program_identification(const std::string& program_name) {
    std::ostringstream s;
    s << program_name;
    s << " (" << get_program_source_identification() << ")";
    return s.str();
}

} // namespace pstrudel

#endif
