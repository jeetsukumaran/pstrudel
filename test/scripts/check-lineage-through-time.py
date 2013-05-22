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
    parser.add_argument("-c", "--check-offsets",
            type=int,
            default=None,
            metavar="NUM-TRANSECTS",
            help="check offsets, assuming NUM-TRANSECTS on tree (default=%(default)s)")
    parser.add_argument("-p", "--precision",
            type=float,
            default=0.01,
            help="numerical precision (default=%(default)s)")
    parser.add_argument("-l", "--label",
            type=str,
            default="[TESTING]",
            help="label for test")
    args = parser.parse_args()

    fails = 0
    stdin_lines = sys.stdin.read().split("\n")
    tree_str = stdin_lines[0]
    tree = dendropy.Tree.get_from_string(tree_str, schema=args.schema)
    # sys.stderr.write("{}: tree read\n".format(args.label))
    reported_offsets = {}
    reported_num_lineages = {}
    for line in stdin_lines[1:]:
        if not line:
            continue
        cols = line.split("\t")
        idx = int(cols[0])
        offset = float(cols[1])
        num_lineages = int(cols[2])
        reported_offsets[idx] = offset
        reported_num_lineages[idx] = num_lineages
    # sys.stderr.write("{}: reported values read\n".format(args.label))
    if not reported_offsets:
        sys.stderr.write("{}: No data on tree reported\n".format(args.label))
        sys.exit(1)
    if args.check_offsets and len(reported_offsets) != args.check_offsets:
        sys.stderr.write("{}: Expecting {} transects but found {}:\n".format(
            args.label,
            args.check_offsets,
            len(reported_offsets),
        ))
        for idx in sorted(reported_offsets.keys()):
            sys.stderr.write("{}, ".format(reported_offsets[idx]))
        sys.stderr.write("\n")
        sys.exit(1)

    max_leaf_distance_from_root = tree.minmax_leaf_distance_from_root()[1]
    offset_step = float(max_leaf_distance_from_root) / args.check_offsets
    for idx in sorted(reported_offsets.keys()):
        offset = reported_offsets[idx]
        calculated_offset = offset_step * (idx + 1)
        if abs(offset - calculated_offset) > args.precision:
            sys.stderr.write("{} Offset distance for index {}: expecting distance of {} but found {}\n".format(
                args.label, idx+1, calculated_offset, offset))
            fails += 1
        reported_num = reported_num_lineages[idx]
        calculated_num = tree.num_lineages_at(offset)
        if reported_num != calculated_num:
            sys.stderr.write("{} Transect at index {}, distance offset from root {}: expecting {} lineages but found {}\n".format(
                args.label, idx+1, offset, calculated_num, reported_num))
            fails += 1
    #     else:
    #         sys.stderr.write("{} Offset at index {}, distance from root {}: OK {} vs. {}\n".format(
    #             args.label, idx+1, offset, calculated_num, reported_num))
    if fails > 0:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == '__main__':
    main()


