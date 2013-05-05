#! /usr/bin/env python

import sys
import dendropy

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
    return subtree_leaf_set_sizes

trees = dendropy.TreeList.get_from_path(sys.argv[1], sys.argv[2])
for tidx, tree in enumerate(trees):
    subtree_leaf_set_sizes = count_subtree_leaf_set_sizes(tree)
    keys = sorted(subtree_leaf_set_sizes.keys())
    for key in keys:
        sys.stdout.write("{}\t{}\t{}\n".format(tidx, key, subtree_leaf_set_sizes[key]))
