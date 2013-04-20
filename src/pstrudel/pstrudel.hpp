#ifndef PSTRUDEL_HPP
#define PSTRUDEL_HPP
#include "version.h"

#include <string>

namespace pstrudel {

std::string get_program_identification(const std::string& program_name) {
    std::ostringstream s;
    s << program_name;
#if defined(PROJECT_GIT_SHA1_SHORT)
    s << " (" << PROJECT_GIT_REFSPEC << ": " << PROJECT_GIT_SHA1_SHORT << ") ";
#else
    s << " (" << __DATE__ << ")";
#endif
    return  s.str();
}

} // namespace pstrudel

#endif
