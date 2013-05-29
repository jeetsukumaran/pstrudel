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
    rep_gamma = float(stdin_lines.pop(0))
    fails += compare_almost_equal(
            args.label,
            "Pybus-Harvey Gamma",
            tree.pybus_harvey_gamma(),
            rep_gamma,
            args.precision)
    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()



