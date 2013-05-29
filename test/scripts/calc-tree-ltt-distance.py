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

def compare_almost_equal(label, comparison_name, expected, observed, check_precision):
    if expected == 0:
        d = 1
    else:
        d = expected
    if abs(observed-expected) > check_precision:
        sys.stderr.write("{}: {}: expecting {:.22f} but found {:.22f}\n".format(label, comparison_name, expected, observed))
        return 1
    else:
        return 0

def compare_almost_equal_vectors(label, comparison_name, expected, observed, check_precision):
    if len(expected) != len(observed):
        sys.stderr.write("{}: {}: expecting profile length of {} but found {}\n".format(label, comparison_name, len(expected), len(observed)))
        for idx in range(max(len(expected), len(observed))):
            if idx >= len(expected):
                ev = "{:24}".format("N/A")
            else:
                ev = "{:.22f}".format(expected[idx])
            if idx >= len(observed):
                ov = "{:24}".format("N/A")
            else:
                ov = "{:.22f}".format(observed[idx])
            sys.stderr.write("    {:03d}: {}    {}\n".format(idx, ev, ov))
        return 1
    fails = 0
    for idx in range(max(len(expected), len(observed))):
        ev = expected[idx]
        ov = observed[idx]
        if abs(ev - ov) > check_precision:
            fails += 1
            leader = "{} FAIL:".format(comparison_name)
            sys.stderr.write("{}   Index {:03d} (Max={:03d}): {:.22f}    {:.22f}\n".format(leader, idx, len(expected)-1, ev, ov))
    return fails

def preprocess_tree(tree, num_transects=None):
    tree_length = tree.length()
    num_leaves = len(tree.leaf_nodes())
    min_edge_length = None
    for edge in tree.postorder_edge_iter():
        if edge.length > 0.0 and (min_edge_length > edge.length or min_edge_length is None):
            min_edge_length = edge.length
    if min_edge_length is None:
        min_edge_length = 0
    max_leaf_distance_from_root = tree.minmax_leaf_distance_from_root()[1]
    max_leaf_distance_from_root = max_leaf_distance_from_root - (float(min_edge_length)/1e3) # account for error
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

def preprocess_trees(trees, num_transects=None):
    for tree in trees:
        preprocess_tree(tree=tree, num_transects=num_transects)

def main():
    parser = argparse.ArgumentParser()
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
            default=1e-14,
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

    stdin_lines = sys.stdin.read().split("\n")
    taxon_set = dendropy.TaxonSet()

    tree1 = dendropy.Tree.get_from_string(stdin_lines.pop(0), args.schema, taxon_set=taxon_set)
    tree2 = dendropy.Tree.get_from_string(stdin_lines.pop(0), args.schema, taxon_set=taxon_set)

    preprocess_tree(tree1)
    preprocess_tree(tree2)

    exp_ltt = euclidean_distance(tree1.lineage_accumulation_through_time, tree2.lineage_accumulation_through_time)
    exp_lst = euclidean_distance(tree1.lineage_splitting_times, tree2.lineage_splitting_times)
    exp_slst = euclidean_distance(tree1.scaled_lineage_splitting_times, tree2.scaled_lineage_splitting_times)

    rep_ltt = float(stdin_lines.pop(0))
    rep_lst = float(stdin_lines.pop(0))
    rep_slst = float(stdin_lines.pop(0))

    rep_tree_profiles = { tree1: {}, tree2: {} }
    rep_tree_profiles[tree1]["ltt"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    rep_tree_profiles[tree1]["lst"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    rep_tree_profiles[tree1]["slst"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    rep_tree_profiles[tree2]["ltt"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    rep_tree_profiles[tree2]["lst"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    rep_tree_profiles[tree2]["slst"] = [float(i) for i in stdin_lines.pop(0).split(",") if i]

    fails = 0
    for tree_idx, tree in enumerate((tree1, tree2)):
        rep_ltt_profile = rep_tree_profiles[tree]["ltt"]
        rep_lst_profile = rep_tree_profiles[tree]["lst"]
        rep_slst_profile = rep_tree_profiles[tree]["slst"]
        fails += compare_almost_equal_vectors(
                args.label,
                "Tree {}: lineage accumulation through time profile".format(tree_idx),
                tree.lineage_accumulation_through_time,
                rep_ltt_profile,
                args.precision)
        fails += compare_almost_equal_vectors(
                args.label,
                "Tree {}: lineage splitting time profile".format(tree_idx),
                tree.lineage_splitting_times,
                rep_lst_profile,
                args.precision)
        fails += compare_almost_equal_vectors(
                args.label,
                "Tree {}: scaled lineage splitting time profile".format(tree_idx),
                tree.scaled_lineage_splitting_times,
                rep_slst_profile,
                args.precision)
    fails += compare_almost_equal(
            args.label,
            "lineage accumulation through time difference",
            exp_ltt,
            rep_ltt,
            args.precision)
    fails += compare_almost_equal(
            args.label,
            "lineage splitting times difference",
            exp_lst,
            rep_lst,
            args.precision)
    fails += compare_almost_equal(
            args.label,
            "scaled lineage splitting times difference",
            exp_slst,
            rep_slst,
            args.precision)

    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()



