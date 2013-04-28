#include <vector>
#include <iostream>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>

int main() {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    pstrudel::treeio::read_from_string(trees, "(a,(b,c)); (c,(b,a));", "newick");
}
