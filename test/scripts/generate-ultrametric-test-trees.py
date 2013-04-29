#! /usr/bin/env python

import os
import sys
import dendropy
import random
import argparse
import itertools
from dendropy import coalescent

"""
Generates ultrametric trees from tree patterns.
"""

def generate_all_permutations(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in generate_all_permutations(elements[1:]):
            for i in range(len(elements)):
                #nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]

def generate_all_labeled_histories(tree,
        pop_size=None,
        rng=None):

    # get collection of nodes organized by the order
    # of coalescence as given by the level
    node_event_order = {}
    event_nodes = {}
    num_leaves = 0
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            node_event_order[nd] = 0
            num_leaves += 1
        else:
            eo = max([node_event_order[ch] for ch in nd.child_nodes()]) + 1
            node_event_order[nd] = eo
        try:
            event_nodes[node_event_order[nd]].append(nd)
        except KeyError:
            event_nodes[node_event_order[nd]] = [nd]
    event_idx_list = sorted(event_nodes.keys())
    event_permutations = []
    for event_idx in event_idx_list[1:]:
        h = list(generate_all_permutations([nd for nd in event_nodes[event_idx]]))
        event_permutations.append(h)
    trees = []
    for idx, ep in enumerate(itertools.product(*event_permutations)):
        # event_nodes_x = { 0: event_nodes[0] }
        # for hidx in range(len(ep)):
        #     event_nodes_x[hidx+1] = ep[hidx]
        coal_event_list = [event_nodes[0]] + list(ep)
        _generate_labeled_history(tree,
                coal_event_list=coal_event_list,
                num_leaves=num_leaves,
                pop_size=pop_size,
                rng=rng)
        labeled_history_tree = dendropy.Tree(tree)
        trees.append(labeled_history_tree)
    return trees

# coal_event_list = a list of lists, with each inner list collecting the nodes
# that all have to coalesce before any node in the next inner list
def _generate_labeled_history(tree,
        coal_event_list,
        num_leaves=None,
        pop_size=None,
        rng=None):
    # set the node ages
    node_ages = {}
    for nd in coal_event_list[0]:
        node_ages[nd] = 0.0
    current_tree_max_age = 0.0
    if num_leaves is None:
        num_nodes_remaining = num_leaves
    else:
        num_nodes_remaining = len([lf for lf in tree.leaf_iter()])
    # process events in order of occurrence
    for event_nodes in coal_event_list[1:]:
        for nd in event_nodes:
            child_nodes = nd.child_nodes()
            max_ch_age = max([node_ages[ch] for ch in child_nodes])
            min_age = max(current_tree_max_age, max_ch_age)
            if rng is None:
                # rng is not given: use expected waiting time
                wt = coalescent.expected_tmrca(
                        n_genes=num_nodes_remaining,
                        pop_size=pop_size,
                        n_to_coalesce=len(child_nodes))
            else:
                # rng is given: draw random waiting time
                wt = coalescent.time_to_coalescence(
                        n_genes=num_nodes_remaining,
                        pop_size=pop_size,
                        n_to_coalesce=len(child_nodes))
            node_ages[nd] = min_age + wt
            current_tree_max_age += wt
            num_nodes_remaining = num_nodes_remaining - len(child_nodes) + 1
    # set the edge lengths according to the node ages
    for nd in tree.preorder_node_iter():
        if nd.parent_node is None:
            nd.edge.length = 0.0
        else:
            nd.edge.length = node_ages[nd.parent_node] - node_ages[nd]
    return tree

def generate_single_labeled_history(tree,
        pop_size=None,
        rng=None):
    # get collection of nodes organized by the order
    # of coalescence as given by the level
    node_event_order = {}
    event_nodes = {}
    num_leaves = 0
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            node_event_order[nd] = 0
            num_leaves += 1
        else:
            eo = max([node_event_order[ch] for ch in nd.child_nodes()]) + 1
            node_event_order[nd] = eo
        try:
            event_nodes[node_event_order[nd]].append(nd)
        except KeyError:
            event_nodes[node_event_order[nd]] = [nd]

    # set the node ages
    node_ages = {}
    event_idx_list = sorted(event_nodes.keys())
    for nd in event_nodes[0]:
        node_ages[nd] = 0.0
    current_tree_max_age = 0.0
    num_nodes_remaining = num_leaves
    for event_idx in event_idx_list[1:]:
        for nd in event_nodes[event_idx]:
            child_nodes = nd.child_nodes()
            max_ch_age = max([node_ages[ch] for ch in child_nodes])
            min_age = max(current_tree_max_age, max_ch_age)
            if rng is None:
                # rng is not given: use expected waiting time
                wt = coalescent.expected_tmrca(
                        n_genes=num_nodes_remaining,
                        pop_size=pop_size,
                        n_to_coalesce=len(child_nodes))
            else:
                # rng is given: draw random waiting time
                wt = coalescent.time_to_coalescence(
                        n_genes=num_nodes_remaining,
                        pop_size=pop_size,
                        n_to_coalesce=len(child_nodes))
            node_ages[nd] = min_age + wt
            current_tree_max_age += wt
            num_nodes_remaining = num_nodes_remaining - len(child_nodes) + 1

    # set the edge lengths according to the node ages
    for nd in tree.preorder_node_iter():
        if nd.parent_node is None:
            nd.edge.length = 0.0
        else:
            nd.edge.length = node_ages[nd.parent_node] - node_ages[nd]
    return tree

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("files", metavar="FILE", type=str, nargs="*", help="path to input data file(s)")
    parser.add_argument("-f", "--input-format",
            default="nexus",
            help="format of input files")
    parser.add_argument("-F", "--output-format",
            default="nexus",
            help="format of output")
    parser.add_argument("-a", "--all-histories",
            action="store_true",
            default=False,
            help="generate *all* labeled histories for each input tree pattern (default is to generate just one)")
    parser.add_argument("-k", "--preserve-labels",
            action="store_true",
            default=False,
            help="preserve original labels (otherwise trees will be relabeled)")
    src_trees = dendropy.TreeList()
    args = parser.parse_args()
    if len(args.files) == 0:
        sys.exit("Need to specify path to one or more sources of trees")
    for fpath in args.files:
        src_trees.read_from_path(fpath, args.input_format)
    result_trees = dendropy.TreeList(taxon_set=src_trees.taxon_set)
    if not args.preserve_labels:
        for tidx, tree in enumerate(src_trees):
            for idx, nd in enumerate(tree.postorder_node_iter()):
                if nd.is_leaf():
                    nd.taxon.label = "t{}".format(idx+1)
                else:
                    nd.label = "i{}".format(idx+1)
    total_tree_count = 0
    for tidx, tree in enumerate(src_trees):
        if args.all_histories:
            gtrees = generate_all_labeled_histories(tree)
            for gidx, gtree in enumerate(gtrees):
                total_tree_count += 1
                gtree.label = "H{} ({}.{})".format(total_tree_count, tidx+1, gidx+1)
            result_trees.extend(gtrees)
        else:
            generate_single_labeled_history(tree)
            total_tree_count += 1
            tree.label = "H{} ({})".format(total_tree_count, tidx+1)
            result_trees.append(tree)
    result_trees.write_to_stream(sys.stdout,
            args.output_format,
            edge_label_compose_func=lambda x: "{:0.2f}".format(x.length))

if __name__ == "__main__":
    main()
