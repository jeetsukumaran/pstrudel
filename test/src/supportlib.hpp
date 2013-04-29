#ifndef PSTRUDEL_TEST_SUPPORT_LIB
#define PSTRUDEL_TEST_SUPPORT_LIB

#include <vector>
#include <iostream>
#include <string>

namespace pstrudel { namespace test {

template <typename S>
void log(S&) {}

template <typename S, typename T>
void log(S& stream, const T& arg1) {
    stream << arg1;
}

template <typename S, typename T>
void log(S& stream, const std::vector<T>& arg1) {
    write_container(arg1, stream, ", ");
    log(stream);
}

template <typename S, typename T, typename... Types>
void log(S& stream, const std::vector<T>& arg1, const Types&... args) {
    write_container(arg1, stream, ", ");
    log(stream, args...);
}

template <typename S, typename T, typename... Types>
void log(S& stream, const T& arg1, const Types&... args) {
    stream << arg1;
    log(stream, args...);
}

template <typename T, typename U, typename... Types>
int fail_test(const std::string& test_name,
        const T& expected,
        const U& observed,
        const Types&... args) {
    std::cerr << "\n||[[ FAIL ]]]";
    std::cerr << "\n||     Test: " << test_name;
    log(std::cerr, "\n|| Expected: ", expected);
    log(std::cerr, "\n|| Observed: ", observed);
    log(std::cerr, "\n||  Remarks: ", args...);
    std::cerr << std::endl;
    return EXIT_FAILURE;
}

} } // namespace pstrudel::test

#endif
