#include <vector>
#include <iostream>
#include <pstrudel/distancetree.hpp>
#include <pstrudel/dataio.hpp>

#include <unistd.h>
#include <stdio.h>
int main() {
    typedef pstrudel::DistanceTree  TreeType;
    std::vector<TreeType>  trees;
    pstrudel::treeio::read_from_filepath(trees, "data/trees/basic-patterns/n06-rooted-patterns.nex", "nexus");
}
