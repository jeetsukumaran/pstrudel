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
    else:
        return 0

def main():
    parser = argparse.ArgumentParser()
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
    tree = dendropy.Tree.get_from_string(stdin_lines.pop(0), args.schema)

    fails = 0

    exp_node_ages = [age for age in tree.calc_node_ages() if age]
    rep_node_ages = [float(i) for i in stdin_lines.pop(0).split(",") if i]
    fails += compare_almost_equal_vectors(
            args.label,
            "Node Ages",
            exp_node_ages,
            rep_node_ages,
            args.precision)

    rep_b1 = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "B1",
            tree.B1(),
            rep_b1,
            args.precision)

    rep_colless = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "Colless' Imbalance",
            tree.colless_tree_imbalance(),
            rep_colless,
            args.precision)

    rep_gamma = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "Pybus-Harvey Gamma",
            tree.pybus_harvey_gamma(),
            rep_gamma,
            args.precision)

    rep_nbar = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "N-Bar",
            tree.N_bar(),
            rep_nbar,
            args.precision)

    rep_treeness = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "treeness",
            tree.treeness(),
            rep_treeness,
            args.precision)

    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()



