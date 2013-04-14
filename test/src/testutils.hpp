#ifndef PSTRUDEL_TESTUTILS_HPP
#define PSTRUDEL_TESTUTILS_HPP

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include "../../src/dataio.hpp"
#include "../../src/pairwise_distance_tree.hpp"
#include "../../src/split.hpp"

template <class TreeT>
void get_trees(std::vector<TreeT>& trees, int argc, const char * argv[]) {
    if (argc == 1) {
        std::cerr << "Usage: " << argv[0] << " <TREEFILE> [nexus|newick]" << std::endl;
        exit(1);
    }
    if (argc > 3) {
        std::cerr << "Expecting at most two arguments" << std::endl;
        exit(1);
    }
    std::string format;
    if (argc == 3) {
        format = argv[2];
    } else {
        format = "nexus";
    }
    std::ifstream src(argv[1]);
    pstrudel::treeio::read_from_stream(trees, src, format);
}

template <typename T>
long read_data_vector(std::istream& is, std::vector<T>& data) {
    std::string tmp;
    if (std::getline(std::cin, tmp).good()) {
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
long read_data_vectors(std::istream& is, std::vector<std::vector<T>>& data) {
    std::string tmp;
    unsigned long count = 0;
    while (std::getline(std::cin, tmp).good()) {
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

#endif
