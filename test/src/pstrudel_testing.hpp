#ifndef PSTRUDEL_TESTTING_HPP
#define PSTRUDEL_TESTTING_HPP

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
// Logging and printing

template<class T>
void write_container(std::ostream& out, const T& container, const std::string& separator=", ") {
    std::copy(container.cbegin(), container.cend(), std::ostream_iterator<typename T::value_type>(out, separator.c_str()));
}

template<class T>
std::string join_container(const T& container, const std::string& separator=", ") {
    std::ostringstream out;
    write_container(out, container, separator);
    return out.str();
}

template <typename S>
void log(S&) {}

template <typename S, typename T>
void log(S& stream, const T& arg1) {
    stream << arg1;
}

template <typename S, typename T>
void log(S& stream, const std::vector<T>& arg1) {
    write_container(stream, arg1, ", ");
    log(stream);
}

template <typename S, typename T, typename... Types>
void log(S& stream, const std::vector<T>& arg1, const Types&... args) {
    write_container(stream, arg1, ", ");
    log(stream, args...);
}

template <typename S, typename T, typename... Types>
void log(S& stream, const T& arg1, const Types&... args) {
    stream << arg1;
    log(stream, args...);
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

//////////////////////////////////////////////////////////////////////////////
// Testing

template <typename T, typename U, typename... Types>
int fail_test(const std::string& test_name,
        unsigned long line_num,
        const T& expected,
        const U& observed,
        const Types&... args) {
    std::cerr << "\n||| FAIL |||";
    std::cerr << "\n|     Test: " << test_name;
    std::cerr << "\n|     Line: " << line_num;
    log(std::cerr, "\n| Expected: ", expected);
    log(std::cerr, "\n| Observed: ", observed);
    log(std::cerr, "\n|  Remarks: ", args...);
    std::cerr << std::endl;
    return EXIT_FAILURE;
}

template <typename T, typename U, typename... Types>
int check_equal(
        const T& expected,
        const U& observed,
        const std::string& test_name,
        unsigned long line_num,
        const Types&... args) {
    if (expected != observed) {
        return fail_test(test_name,
                line_num,
                expected,
                observed,
                args...);
    } else {
        return EXIT_SUCCESS;
    }
}

inline bool is_almost_equal(double a, double b, double tolerance=1e-6) {
//     if (std::abs(a-b)<2.*std::numeric_limits<double>::epsilon() || // For small a, b
//             std::abs(a-b)<std::min(a,b)/1.e15) {
//         return true;
//     } else {
//         return false;
//     }
    double abs_a = std::fabs(a);
    double abs_b = std::fabs(b);
    double diff = std::fabs(a - b);
    const static double MIN_DOUBLE_VALUE = std::numeric_limits<double>::min();
    if (a == b) {
        // shortcut, handles infinities
        return true;
    } else if (abs_a < tolerance && abs_b < tolerance) {
        return true;
    } else if (a == 0 || b == 0 || diff < MIN_DOUBLE_VALUE)  {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        return diff < (tolerance * MIN_DOUBLE_VALUE);
    } else {
        // use relative error
        return (diff / (abs_a + abs_b)) < tolerance;
    }

}

template <typename... Types>
int check_almost_equal(
        double expected,
        double observed,
        const std::string& test_name,
        unsigned long line_num,
        const Types&... args) {
    if (!is_almost_equal(expected, observed))  {
        return fail_test(test_name,
                line_num,
                expected,
                observed,
                args...);
    } else {
        return EXIT_SUCCESS;
    }
}

} } // namespace pstrudel::test

#endif
