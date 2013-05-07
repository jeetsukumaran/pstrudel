#! /usr/bin/env python

import sys
import math
import collections
import dendropy

def count_set_size_difference(v1, v2):
    c1 = collections.Counter(v1)
    c2 = collections.Counter(v2)
    counted_matched = c1 & c2
    matched = sorted(list(counted_matched.elements()))
    counted_diffs = (c1 - c2) + (c2 - c1)
    unmatched = sorted(list(counted_diffs.elements()))
    diff = len(unmatched)
    return diff

def count_subtree_leaf_set_sizes(tree):
    internal_nodes = tree.internal_nodes()
    subtree_leaf_set_sizes = {}
    for nd in internal_nodes:
        leaf_count = 0
        for leaf in nd.leaf_iter():
            leaf_count += 1
        if nd.taxon is not None:
            label = nd.taxon.label
        else:
            label = nd.label
        subtree_leaf_set_sizes[label] = leaf_count
    return sorted(subtree_leaf_set_sizes.values())

def main():
    trees = dendropy.TreeList.get_from_path(sys.argv[1], sys.argv[2])
    for tree in trees:
        tree.subtree_leaf_set_sizes = count_subtree_leaf_set_sizes(tree)
    for tidx1 in range(len(trees)):
        for tidx2 in range(len(trees)):
            sys.stdout.write("{}\t{}\t{}\n".format(
                tidx1,
                tidx2,
                count_set_size_difference(
                        trees[tidx1].subtree_leaf_set_sizes,
                        trees[tidx2].subtree_leaf_set_sizes)))

if __name__ == "__main__":
    main()



