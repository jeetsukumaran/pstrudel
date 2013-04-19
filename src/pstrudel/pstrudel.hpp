#ifndef PSTRUDEL_HPP
#define PSTRUDEL_HPP

#include <string>

namespace pstrudel {

std::string get_program_identification(const std::string& program_name) {
    std::ostringstream s;
    s << program_name;
    // s << " (" << PACKAGE_NAME << " v" << PACKAGE_VERSION << ": ";
#if defined(PSTRUDEL_SOURCE_DESC)
    s << PSTRUDEL_SOURCE_DESC << ")";
#else
    s << __DATE__ << ")";
#endif
    return  s.str();
}

} // namespace pstrudel

#endif