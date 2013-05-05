#ifndef PSTRUDEL_TESTTING_HPP
#define PSTRUDEL_TESTTING_HPP

#include <platypus/utility/testing.hpp>
#include <iomanip>
#include <limits>
#include <vector>
#include <sstream>
#include <iterator>
#include <iostream>
#include <string>
#include <cmath>

#ifdef _WIN32
const char FILESYSTEM_PATH_SEP_CHAR = '\\';
#else
const char FILESYSTEM_PATH_SEP_CHAR = '/';
#endif
const std::string FILESYSTEM_PATH_SEP_STR = std::string(1, FILESYSTEM_PATH_SEP_CHAR);

namespace pstrudel { namespace test {

//////////////////////////////////////////////////////////////////////////////
// Reading

template <typename T>
long read_data_vector(std::istream& input, std::vector<T>& data) {
    std::string tmp;
    if (std::getline(input, tmp).good()) {
        std::stringstream ss(tmp);
        unsigned long count = 0;
        T val;
        while (ss >> val) {
            data.push_back(val);
            ++count;
        }
        return count;
    } else {
        return -1;
    }
}

template <typename T>
long read_data_vectors(std::istream& input, std::vector<std::vector<T>>& data) {
    std::string tmp;
    unsigned long count = 0;
    while (std::getline(input, tmp).good()) {
        data.emplace_back();
        auto & v = data.back();
        std::stringstream ss(tmp);
        T val;
        while (ss >> val) {
            v.push_back(val);
        }
        ++count;
    }
    return count;
}

//////////////////////////////////////////////////////////////////////////////
// Paths and filesystem

std::string get_dir(const char * argv0);
std::string get_test_dir(const char * argv0);

template <typename T, typename S>
std::string join_path(const T & arg1, const S & arg2) {
    return arg1 + FILESYSTEM_PATH_SEP_STR + arg2;
}

template <typename T, typename ... Types>
std::string join_path(const T & arg1, const Types & ... args) {
    return arg1 + FILESYSTEM_PATH_SEP_STR + join_path(args ...);
}

//////////////////////////////////////////////////////////////////////////////
// External process execution

std::string execute_external_process(const std::string& cmd,
        bool error_on_error_exit=true,
        bool error_on_timeout_empty=true);
        // unsigned long timeout=1000);

} } // namespace pstrudel::test

#endif
