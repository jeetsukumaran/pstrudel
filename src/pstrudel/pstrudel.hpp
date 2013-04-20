#ifndef PSTRUDEL_HPP
#define PSTRUDEL_HPP
#include "version.h"

#include <string>

namespace pstrudel {

std::string get_program_identification(const std::string& program_name) {
    std::ostringstream s;
    s << program_name;
#if defined(PROJECT_SOURCE_DESC)
    s << " (" PROJECT_SOURCE_DESC << ") ";
#else
    s << " (" << __DATE__ << ")";
#endif
    return  s.str();
}

} // namespace pstrudel

#endif
