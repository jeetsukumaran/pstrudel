#! /usr/bin/env python

###############################################################################
##
##  Copyright 2012 Jeet Sukumaran.
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License along
##  with this program. If not, see <http://www.gnu.org/licenses/>.
##
###############################################################################

import sys
import os
import argparse
import dendropy
from dendropy import coalescent

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
    parser.add_argument("-l", "--label",
            type=str,
            default="[TESTING]",
            help="label for test")
    parser.add_argument("-r", "--regime",
            dest="regime",
            type=str,
            default="mean",
            choices=["mean", "uniform", "anti-coalescent"],
            help="expected regime under which the coalescent intervals were generated (default='%(default)s')")
    parser.add_argument("-p", "--precision",
            type=float,
            default=0.01,
            help="numerical precision (default=%(default)s)")
    args = parser.parse_args()

    fails = 0
    trees = dendropy.TreeList.get_from_stream(sys.stdin, schema=args.schema)
    for tidx, tree in enumerate(trees):

        # check structure
        if not tree._debug_tree_is_valid():
            sys.stderr.write("{} Tree {}: structure is invalid:\n{}".format(args.label, tidx+1, tree.compose_newick()))
            fails += 1

        # ensure parent ages > child ages
        tree.calc_node_ages(check_prec=args.precision)
        for nd in tree.postorder_node_iter():
            if nd.parent_node is not None:
                if nd.age > nd.parent_node.age:
                    sys.stderr.write("{} Tree {}: Node '{}': age ({}) is greater than parent age ({})".format(
                        args.label, tidx+1, nd.label, nd.age, nd.parent_node.age))
                    fails += 1

        # check waiting times
        num_tips = len(tree.leaf_nodes())
        wf = coalescent.extract_coalescent_frames(tree, check_ultrametricity_prec=args.precision)
        num_lineages = sorted(wf.keys())
        for n in num_lineages:
            wt = wf[n]
            if args.regime == "mean":
                exp_wt = coalescent.expected_tmrca(n)
            elif args.regime == "uniform":
                exp_wt = 1.0
            elif args.regime == "anti-coalescent":
                exp_wt = coalescent.expected_tmrca(num_tips - n)
            if abs(exp_wt - wt) > args.precision:
                sys.stderr.write("{} Tree {}: Waiting time for coalescence event with {} lineages: expecting {} but found {}\n".format(
                        args.label, tidx+1, n, exp_wt, wt))
                fails += 1
    sys.exit(1)
    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == '__main__':
    main()


