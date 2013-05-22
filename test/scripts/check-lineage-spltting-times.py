#! /usr/bin/env python

import sys
import os
import argparse
import dendropy

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--format",
            dest="schema",
            type=str,
            default="nexus",
            choices=["nexus", "newick"],
            help="input data format (default='%(default)s')")
    parser.add_argument("-p", "--precision",
            type=float,
            default=1e-8,
            help="numerical precision (default=%(default)s)")
    parser.add_argument("-l", "--label",
            type=str,
            default="[TESTING]",
            help="label for test")
    args = parser.parse_args()

    fails = 0
    message_prefix = "[{}]".format(args.label)
    stdin_lines = sys.stdin.read().split("\n")
    tree_str = stdin_lines[0]
    tree = dendropy.Tree.get_from_string(tree_str, schema=args.schema)
    reported_splitting_times = []
    for line in stdin_lines[1:]:
        if not line:
            continue
        st = float(line)
        reported_splitting_times.append(st)
    tree.calc_node_root_distances()
    expected_splitting_times = []
    for node in tree.preorder_node_iter():
        if not node.is_leaf():
            expected_splitting_times.append(node.root_distance)
    expected_splitting_times.sort()
    if len(expected_splitting_times) != len(reported_splitting_times):
            sys.stderr.write("{} Expected {} splitting times but found {}\n".format(
                message_prefix, len(expected_splitting_times), len(reported_splitting_times)))
    for idx, exp_st in enumerate(expected_splitting_times):
        rep_st = reported_splitting_times[idx]
        if abs(rep_st - exp_st) > args.precision:
            sys.stderr.write("{} Time for split event {}: expecting {} but found {}\n".format(
                message_prefix, idx+1, exp_st, rep_st))
            fails += 1
    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == '__main__':
    main()


