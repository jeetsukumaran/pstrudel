#include <sstream>
#include <stdexcept>
#include "supportlib.hpp"

namespace pstrudeltest {
const unsigned long MAX_PROCESS_BUFFER_SIZE = 10000;

std::string dir_path(const char * argv0) {
    std::string dpath = argv0;
    dpath = dpath.substr(0, dpath.find_last_of(FILESYSTEM_PATH_SEP));
    return dpath;
}

std::string execute_external_process(const std::string& cmd,
        bool error_on_error_exit,
        bool error_on_timeout_empty,
        unsigned long timeout) {
    static char stdout_buffer[MAX_PROCESS_BUFFER_SIZE];
    // LOG.debug() << "Executing: '" << cmd << "'" << std::endl;
    FILE * process_stdout = popen(cmd.c_str(), "r");
    unsigned long max_wait_cycles = timeout;
    stdout_buffer[0] = '\0';
    while (stdout_buffer[0] == '\0' && max_wait_cycles > 0) {
        fgets(stdout_buffer, 1000, process_stdout);
        --max_wait_cycles;
    }
    int ret_code = pclose(process_stdout);
    if (stdout_buffer[0] == '\0' && error_on_timeout_empty) {
        std::cerr << "\nFailed to read data from standard output of process:\n\n   ";
        std::cerr << cmd << "\n" << std::endl;
        throw std::runtime_error("External process failed");
    } else if (ret_code != 0) {
        std::cerr << "\nProcess returned non-zero exit code " << ret_code << ":\n\n   ";
        std::cerr << cmd << "\n" << std::endl;
        std::cerr << stdout_buffer << std::endl;
        throw std::runtime_error("External process failed");
    }
    return std::string(stdout_buffer);
}




} // namespace pstrudeltest
