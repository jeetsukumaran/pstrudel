#! /usr/bin/env python

import os
import sys
import math
import argparse
import dendropy
from dendropy import treecalc

def euclidean_distance(v1, v2):
    ss = 0.0
    assert len(v1) == len(v2)
    for idx1 in range(len(v1)):
        ss += pow((v1[idx1] - v2[idx1]), 2)
    return math.sqrt(ss)

def preprocess_trees(trees, num_transects=None):
    for tree in trees:
        tree_length = tree.length()
        num_leaves = len(tree.leaf_nodes())
        max_leaf_distance_from_root = tree.minmax_leaf_distance_from_root()[1]
        if not num_transects:
            num_transects = (num_leaves - 1) * 10
        offset_step = float(max_leaf_distance_from_root) / num_transects
        transect_offsets = []
        offset = offset_step
        for transect_idx in range(num_transects):
            transect_offsets.append(offset)
            offset += offset_step
        tree.lineage_accumulation_through_time = []
        for tof in transect_offsets:
            if tof > max_leaf_distance_from_root:
                tof = max_leaf_distance_from_root
            n = float(tree.num_lineages_at(tof))/num_leaves
            tree.lineage_accumulation_through_time.append(n)
        tree.lineage_splitting_times = []
        tree.scaled_lineage_splitting_times = []
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                continue
            if tree_length > 0:
                tree.lineage_splitting_times.append(node.root_distance)
                tree.scaled_lineage_splitting_times.append(float(node.root_distance)/tree_length)
            else:
                tree.lineage_splitting_times.append(0.0)
                tree.scaled_lineage_splitting_times.append(0.0)
        tree.lineage_splitting_times.sort()
        tree.scaled_lineage_splitting_times.sort()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("treefile",
            type=str,
            nargs="*",
            help="path to tree file(s)")
    parser.add_argument("-",
            action="store_true",
            default=False,
            dest="trees_from_stdin",
            help="read trees from standard input")
    parser.add_argument("-f", "--format",
            dest="schema",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="input data format (default='%(default)s')")
    parser.add_argument("-n", "--num-transects",
            type=int,
            default=None,
            help="number of transects or slices on tree for lineage accumulation counts (default=%(default)s)")
    parser.add_argument("-p", "--precision",
            type=int,
            default=22,
            help="numerical precision (default=%(default)s)")
    parser.add_argument("-l", "--label",
            type=str,
            default="[TESTING]",
            help="label for test")
    parser.add_argument("--verbosity",
            type=int,
            default=0,
            help="level of messaging")
    args = parser.parse_args()

    trees = dendropy.TreeList()
    if args.trees_from_stdin:
        trees.read_from_stream(sys.stdin, schema=args.schema)
    for fpath in args.treefile:
        trees.read_from_path(os.path.expanduser(os.path.expandvars(fpath)), schema=args.schema)
    if len(trees) == 0:
        sys.exit("Need to specify '-' option and pipe trees to standard input and/or specify paths to at least one tree file as argument")
    preprocess_trees(trees)
    out_template = "{{}}\t{{}}\t{{:.{precision}f}}\t{{:.{precision}f}}\t{{:.{precision}f}}\n".format(precision=args.precision)
    # out_template = "{}\t{}\t{:f}\t{:f}\t{:f}\n"
    for tidx1 in range(len(trees)):
        if (args.verbosity > 1):
            sys.stderr.write("Tree {}: {}\n".format(tidx1, ", ".join("{:.22f}".format(i) for i in trees[tidx1].lineage_accumulation_through_time)))
            sys.stderr.write("Tree {}: {}\n".format(tidx1, ", ".join("{:.22f}".format(i) for i in trees[tidx1].lineage_splitting_times)))
            sys.stderr.write("Tree {}: {}\n".format(tidx1, ", ".join("{:.22f}".format(i) for i in trees[tidx1].scaled_lineage_splitting_times)))
        for tidx2 in range(len(trees)):
            sys.stdout.write(out_template.format(
                tidx1,
                tidx2,
                euclidean_distance(trees[tidx1].lineage_accumulation_through_time , trees[tidx2].lineage_accumulation_through_time) ,
                euclidean_distance(trees[tidx1].lineage_splitting_times           , trees[tidx2].lineage_splitting_times          )  ,
                euclidean_distance(trees[tidx1].scaled_lineage_splitting_times    , trees[tidx2].scaled_lineage_splitting_times   )
                ))


if __name__ == "__main__":
    main()



