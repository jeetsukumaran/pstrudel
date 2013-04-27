#! /usr/bin/env python

"""
Run tests.
"""

import sys
import os
import argparse
import inspect
import subprocess
import random
import math
import collections
import re
from cStringIO import StringIO

import dendropy
from dendropy import treecalc
from dendropy.interop import paup

class AnsiColorMeta(type):

    ##############################################################################
    ## Color infrastructure modified from:
    ##
    ## http://dev.pocoo.org/hg/pygments-main/file/b2deea5b5030/pygments/console.py
    ##
    ## pygments.console
    ## ~~~~~~~~~~~~~~~~
    ## Format colored console output.
    ## :copyright: Copyright 2006-2009 by the Pygments team, see AUTHORS.
    ## :license: BSD, see LICENSE for details.
    ##

    ansiesc = "\x1b["

    @staticmethod
    def get_ansicodes():
        ansicodes = {}
        ansicodes[""]          = ""
        ansicodes["reset"]     = AnsiColorMeta.ansiesc + "39;49;00m"
        ansicodes["bold"]      = AnsiColorMeta.ansiesc + "01m"
        ansicodes["faint"]     = AnsiColorMeta.ansiesc + "02m"
        ansicodes["standout"]  = AnsiColorMeta.ansiesc + "03m"
        ansicodes["underline"] = AnsiColorMeta.ansiesc + "04m"
        ansicodes["blink"]     = AnsiColorMeta.ansiesc + "05m"
        ansicodes["overline"]  = AnsiColorMeta.ansiesc + "06m"
        dark_colors  = ["black", "darkred", "darkgreen", "brown", "darkblue",
                        "purple", "teal", "lightgray"]
        light_colors = ["darkgray", "red", "green", "yellow", "blue",
                        "fuchsia", "turquoise", "white"]
        x = 30
        for d, l in zip(dark_colors, light_colors):
            ansicodes[d] = AnsiColorMeta.ansiesc + "%im" % x
            ansicodes[l] = AnsiColorMeta.ansiesc + "%i;01m" % x
            x += 1
        ansicodes["darkteal"]   = ansicodes["turquoise"]
        ansicodes["darkyellow"] = ansicodes["brown"]
        ansicodes["fuscia"]     = ansicodes["fuchsia"]
        # ansicodes["white"]      = ansicodes["bold"]
        return ansicodes

    def reset_color(cls):
        return cls.ansicodes["reset"]

    def colorize(cls, color_key, text):
        return cls.ansicodes[color_key] + text + cls.ansicodes["reset"]

    def ansiformat(cls, attr, text):
        """
        Format ``text`` with a color and/or some attributes::

            color       normal color
            *color*     bold color
            _color_     underlined color
            +color+     blinking color
        """
        result = []
        if attr[:1] == attr[-1:] == '+':
            result.append(cls.ansicodes['blink'])
            attr = attr[1:-1]
        if attr[:1] == attr[-1:] == '*':
            result.append(cls.ansicodes['bold'])
            attr = attr[1:-1]
        if attr[:1] == attr[-1:] == '_':
            result.append(cls.ansicodes['underline'])
            attr = attr[1:-1]
        result.append(cls.ansicodes[attr])
        result.append(text)
        result.append(cls.ansicodes['reset'])
        return ''.join(result)

    def __new__(cls, name, bases, dct):
        return type.__new__(cls, name, bases, dct)

    def __init__(cls, name, bases, dct):
        super(AnsiColorMeta, cls).__init__(name, bases, dct)
        # setattr(cls, "ansicodes", AnsiColorMeta.get_ansicodes())
        # setattr(cls, "ansiesc", AnsiColorMeta.ansiesc)
        cls.ansicodes = AnsiColorMeta.get_ansicodes()
        cls.ansiesc = AnsiColorMeta.ansiesc

class AnsiColor(object):

    __metaclass__ = AnsiColorMeta

    def __init__(self, stream=sys.stdout, colorize=True):
        self.stream = stream
        self.colorize = colorize
        self.color_pattern = re.compile(r"@(\w+)@<<(.*?)>>")

    def format_color(self, message):
        if not self.colorize or not self.color_pattern.findall(message):
            return message
        else:
            output = []
            cur_pos = 0
            for match in self.color_pattern.finditer(message):
                start, end = match.span()
                output.append(message[cur_pos:start])
                output.append(AnsiColor.ansiformat(match.group(1), match.group(2)))
                cur_pos = end
            output.append(message[cur_pos:])
            output = "".join(output)
            return output

    def write(self, message):
        self.stream.write(self.format_color(message))

    def __call__(self, message):
        self.write(message)

class TestRunner(object):

    PASS = 0
    FAIL = 1
    ERROR = 2

    dna_to_bases =  {
        'A' :  0,
        'C' :  1,
        'G' :  2,
        'T' :  3,
        'U' :  3,
        'N' :  4,
        'X' :  4,
        '-' :  4,
        '?' :  4,
        'R' :  5,
        'Y' :  6,
        'M' :  7,
        'W' :  8,
        'S' :  9,
        'K' : 10,
        'V' : 11,
        'H' : 12,
        'D' : 13,
        'B' : 14
    }

    dna_to_partials =  {
        'A' : [1.0, 0.0, 0.0, 0.0],
        'C' : [0.0, 1.0, 0.0, 0.0],
        'G' : [0.0, 0.0, 1.0, 0.0],
        'T' : [0.0, 0.0, 0.0, 1.0],
        'U' : [0.0, 0.0, 0.0, 1.0],
        'N' : [1.0, 1.0, 1.0, 1.0],
        'X' : [1.0, 1.0, 1.0, 1.0],
        '-' : [1.0, 1.0, 1.0, 1.0],
        '?' : [1.0, 1.0, 1.0, 1.0],
        'R' : [1.0, 0.0, 1.0, 0.0],
        'Y' : [0.0, 1.0, 0.0, 1.0],
        'M' : [1.0, 1.0, 0.0, 0.0],
        'W' : [1.0, 0.0, 0.0, 1.0],
        'S' : [0.0, 1.0, 1.0, 0.0],
        'K' : [0.0, 0.0, 1.0, 1.0],
        'V' : [1.0, 1.0, 1.0, 0.0],
        'H' : [1.0, 1.0, 0.0, 1.0],
        'D' : [1.0, 0.0, 1.0, 1.0],
        'B' : [0.0, 1.0, 1.0, 1.0]
    }

    @staticmethod
    def get_node_tag(nd):
        if nd.taxon:
            return nd.taxon.label
        else:
            return nd.label

    def __init__(self, opts):
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        self.data_dir = os.path.join(self.script_dir, "data")
        self.verbosity = opts.verbosity
        self.break_on_fail = opts.break_on_fail
        self.random_seed = opts.random_seed
        self.rng = random.Random(self.random_seed)
        self.test_command = None
        self.test_retcode = None
        self.test_stdout = None
        self.test_stderr = None
        self.test_result = None
        self.test_fail_message = None
        self.test_pass_message = None
        self.cout = AnsiColor(sys.stdout, colorize=True)

    def execute_test(self,
            test_program,
            args=None,
            stdin_src=None):
        self.test_command = None
        self.test_stdout = None
        self.test_stderr = None
        self.test_retcode = None
        self.test_fail_message = None
        cmd = [os.path.abspath(os.path.join(self.script_dir, test_program))]
        if args:
            cmd.extend(str(c) for c in args)
        self.test_command = " ".join([str(c) for c in cmd])
        if stdin_src is not None:
            stdin = subprocess.PIPE
            if not stdin_src.endswith("\n"):
                stdin_src += "\n"
        else:
            stdin = None
        p = subprocess.Popen(cmd,
                stdin=stdin,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                cwd=self.script_dir)
        self.test_stdout, self.test_stderr = p.communicate(stdin_src)
        self.test_retcode = p.returncode
        return self.test_retcode

    def fail(self, message):
        self.test_fail_message = message
        return TestRunner.FAIL

    def is_almost_equal(self, v1, v2, prec=1e-4):
        if v1 is None and v2 is None:
            return True
        if v1 is None:
            return self.is_almost_equal(v2, 0.0)
        if v2 is None:
            return self.is_almost_equal(v1, 0.0)
        if abs(v1-v2) > prec:
            return False
        return True

    def compare_trees(self, tree1, tree2):
        status = TestRunner.PASS
        tree1.update_splits()
        tree2.update_splits()

        splits = set(tree1.split_edges.keys() + tree2.split_edges.keys())
        for split in splits:
            if split not in tree2.split_edges:
                return self.fail("Split {} not found on tree 2: {}".format(split, tree1.taxon_set.split_as_newick_string(split)))
            if split not in tree1.split_edges:
                return self.fail("Split {} not found on tree 1: {}".format(split, tree2.taxon_set.split_as_newick_string(split)))
            edge1_len = tree1.split_edges[split].length
            edge2_len = tree2.split_edges[split].length
            if not self.is_almost_equal(edge1_len, edge2_len):
                return self.fail("Unequal edge length for split {} {}: {} vs. {}".format(
                    split,
                    tree1.taxon_set.split_as_newick_string(split),
                    edge1_len,
                    edge2_len))
        if len(tree1.split_edges) != len(tree2.split_edges):
            return self.fail("Different number of splits on trees: {} vs. {}".format(len(tree1.split_edges), len(tree2.split_edges)))
        return status

    def compare_tree_traversal(self, tree1, tree2, traverse_func):
        status = TestRunner.PASS
        tree1_nodes = [nd for nd in getattr(tree1, traverse_func)()]
        tree2_nodes = [nd for nd in getattr(tree2, traverse_func)()]
        if len(tree1_nodes) != len(tree2_nodes):
            return self.fail("Trees have different number of nodes: {} vs. {}".format(len(tree1_nodes), len(tree2_nodes)))
        for nd_idx, node1 in enumerate(tree1_nodes):
            node2 = tree2_nodes[nd_idx]
            if node1.taxon is not node2.taxon:
                return self.fail("Different taxa found during postorder traversal of nodes: {} vs. {}".format(node1.taxon, node2.taxon))
            if node1.label is not node2.label:
                return self.fail("Different labels found during postorder traversal of nodes: {} vs. {}".format(node1.label, node2.label))
            if not self.is_almost_equal(node1.edge.length, node2.edge.length):
                return self.fail("Different edge lengths found during postorder traversal of nodes: {} vs. {}".format(node1.edge.length, node2.edge.length))
        return status

    def test_tree_postorder_iter(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
        self.execute_test("tree_postorder_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_visits = []
        test_edge_lens = []
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            node, edge_len = item.split("\t")
            test_visits.append(node)
            test_edge_lens.append(float(edge_len))
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_visits = []
        check_edge_lens = []
        for check_node in check_tree.postorder_node_iter():
            label = self.get_node_tag(check_node)
            edge_len = check_node.edge.length if check_node.edge.length else 0.0
            check_visits.append(label)
            check_edge_lens.append(edge_len)
        for idx in range(len(check_visits)):
            if idx > len(test_visits):
                return self.fail("Insufficient visits: expecting {} but found {}".format(len(check_visits), len(test_visits)))
            n1 = check_visits[idx]
            n2 = test_visits[idx]
            if n1 != n2:
                return self.fail("Incorrect visit {}: '{}' vs. '{}'".format(idx+1, n1, n2))
            e1 = check_edge_lens[idx]
            e2 = test_edge_lens[idx]
            if not self.is_almost_equal(e1, e2):
                return self.fail("Incorrect node edge length on visit {}: {} vs. {}".format(idx+1, e1, e2))
        return TestRunner.PASS

    def test_tree_leaf_iter(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
        self.execute_test("tree_leaf_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_leaves = []
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            idx, label = item.split("\t")
            test_leaves.append(label)
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_leaves = [ self.get_node_tag(nd) for nd in check_tree.leaf_iter() ]
        test_leaves.sort()
        check_leaves.sort()
        if set(check_leaves) != set(test_leaves):
            return self.fail("Unequal leaf set: {} vs. {}".format(set(check_leaves), set(test_leaves)))
        if len(check_leaves) != len(test_leaves):
            return self.fail("Duplicate leaves: {}".format(set(test_leaves)))
        return TestRunner.PASS

    def test_tree_child_iter(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
        self.execute_test("tree_child_iter",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        test_node_children = {}
        for item in self.test_stdout.split("\n"):
            if not item:
                continue
            nodes = item.split("\t")
            test_node_children[nodes[0]] = nodes[1:]
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", taxon_set=taxa)
        check_node_children = {}
        for check_node in check_tree.postorder_node_iter():
            check_node_children[self.get_node_tag(check_node)] = [self.get_node_tag(child) for child in check_node.child_nodes()]
        for check_node, check_children in check_node_children.items():
            if check_node not in test_node_children:
                return self.fail("Node not visited: '{}'".format(check_node))
            if test_node_children[check_node] != check_children:
                return self.fail("Incorrect children: '{}' vs. '{}'".format(check_children, test_node_children[check_node]))
        return TestRunner.PASS

    def test_tree_reader(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.reference-trees.nexus")
        taxa = dendropy.TaxonSet()
        check_trees = dendropy.TreeList.get_from_path(treefile, "nexus", taxon_set=taxa)
        for code_branch in range(1):
            self.execute_test("tree_reader", [treefile, code_branch, "nexus"])
            if self.test_retcode != 0:
                return TestRunner.ERROR
            test_trees = dendropy.TreeList.get_from_string(self.test_stdout, "newick", taxon_set=taxa)
            if len(check_trees) != len(test_trees):
                return self.fail("TreeReader branch {}: Expecting {} trees in output but found {}".format(
                        code_branch,
                        len(check_trees),
                        len(test_trees)))
                for tidx, tree1 in enumerate(check_trees):
                    status = self.compare_trees(tree1, test_trees[tidx])
                    if status != TestRunner.PASS:
                        return status
        return TestRunner.PASS

    def test_tree_read_from_file(self):
        treefile = os.path.join(self.data_dir, "general", "bird_orders.nex")
        self.execute_test("read_tree",
                [treefile, "nexus"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        taxa = dendropy.TaxonSet()
        check_tree = dendropy.Tree.get_from_path(treefile, "nexus", taxon_set=taxa)
        test_tree = dendropy.Tree.get_from_string(self.test_stdout, "newick", taxon_set=taxa)
        # return self.compare_tree_traversal(check_tree, test_tree, "postorder_node_iter")
        return self.compare_trees(check_tree, test_tree)

    def test_tree_pairwise_tip_distance1(self):
        treefile = os.path.join(self.data_dir, "general", "smdist.tre")
        self.execute_test("tree_pairwise_tip_distances",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        rows = self.test_stdout.split("\n")
        expected = {
                frozenset(['A' ,  'A']) : (  0 , 0  ),
                frozenset(['A' ,  'B']) : (  2 , 3  ),
                frozenset(['A' ,  'C']) : (  5 , 27 ),
                frozenset(['A' ,  'D']) : (  4 , 19 ),
                frozenset(['A' ,  'E']) : (  5 , 26 ),
                frozenset(['B' ,  'B']) : (  0 , 0  ),
                frozenset(['B' ,  'C']) : (  5 , 28 ),
                frozenset(['B' ,  'D']) : (  4 , 20 ),
                frozenset(['B' ,  'E']) : (  5 , 27 ),
                frozenset(['C' ,  'C']) : (  0 , 0  ),
                frozenset(['C' ,  'D']) : (  3 , 14 ),
                frozenset(['C' ,  'E']) : (  2 , 9  ),
                frozenset(['D' ,  'D']) : (  0 , 0  ),
                frozenset(['D' ,  'E']) : (  3 , 13 ),
                frozenset(['E' ,  'E']) : (  0 , 0  ),
                }
        for row in rows:
            if not row:
                continue
            label1, label2, obs_steps, obs_weight = row.split("\t")
            obs_weight = float(obs_weight)
            obs_steps = int(obs_steps)
            exp_steps, exp_weight = expected[frozenset([label1, label2])]
            if not self.is_almost_equal(obs_weight, exp_weight):
                return self.fail("Different path weight for '{}' to '{}': expecting {} but found {}".format(
                    label1,
                    label2,
                    exp_weight,
                    obs_weight))
            if exp_steps != obs_steps:
                return self.fail("Different path step count for '{}' to '{}': expecting {} but found {}".format(
                    label1,
                    label2,
                    exp_steps,
                    obs_steps))
        return TestRunner.PASS

    def test_tree_pairwise_tip_distance2(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
        self.execute_test("tree_pairwise_tip_distances",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        check_tree = dendropy.Tree.get_from_path(treefile, "newick", as_rooted=True)
        pdm = treecalc.PatristicDistanceMatrix(check_tree)
        rows = self.test_stdout.split("\n")
        for row in rows:
            if not row:
                continue
            label1, label2, obs_steps, obs_weight = row.split("\t")
            obs_steps = int(obs_steps)
            obs_weight = float(obs_weight)
            taxon1 = check_tree.taxon_set.get_taxon(label=label1)
            taxon2 = check_tree.taxon_set.get_taxon(label=label2)
            exp_steps = pdm.path_edge_count(taxon1, taxon2)
            exp_weight = pdm(taxon1, taxon2)
            if exp_steps != obs_steps:
                return self.fail("Different path step count for '{}' to '{}': expecting {} but found {}".format(
                    label1,
                    label2,
                    exp_steps,
                    obs_steps))
            if not self.is_almost_equal(obs_weight, exp_weight):
                return self.fail("Different path weight for '{}' to '{}': expecting {} but found {}".format(
                    label1,
                    label2,
                    exp_weight,
                    obs_weight))
            # if exp_steps != obs_steps:
            #     print("Different path step count for '{}' to '{}': expecting {} but found {}".format(
            #         label1,
            #         label2,
            #         exp_steps,
            #         obs_steps))
            # if not self.is_almost_equal(obs_weight, exp_weight):
            #     print("Different path weight for '{}' to '{}': expecting {} but found {}".format(
            #         label1,
            #         label2,
            #         exp_weight,
            #         obs_weight))
        return TestRunner.PASS

    def test_read_dna_sequences(self):
        datafile = os.path.join(self.data_dir, "general", "pythonidae.chars.fasta")
        self.execute_test("read_dna_sequences",
                [datafile, "fasta"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        dna_matrix = dendropy.DnaCharacterMatrix.get_from_path(datafile, "fasta", row_type="STR")
        expected_bases = {}
        for taxon in dna_matrix:
            # label = taxon.label.replace(" ", "_")
            label = taxon.label
            expected_bases[label] = [self.dna_to_bases[state] for state in dna_matrix[taxon]]
            # states = dna_matrix[taxon]
            # for state in states:
            #     sub_bases = self.dna_to_bases[state]
            #     expected_bases[label].extend(list(sub_bases))
        rows = self.test_stdout.split("\n")
        observed_bases = {}
        for row in rows:
            if not row:
                continue
            label, bases = row.split(":")
            bases = [int(v) for v in bases.split(";") if v]
            observed_bases[label] = list(bases)
        for label in expected_bases:
            if label not in observed_bases:
                return self.fail("Sequence '{}' not found: {}".format(label,
                    ",".join(["'{}'".format(t) for t in observed_bases])))
            p1 = expected_bases[label]
            p2 = observed_bases[label]
            if len(p1) != len(p2):
                return self.fail("Sequence '{}': expected {} elements but found {}".format(label,
                    len(p1), len(p2)))
            for idx, i1 in enumerate(p1):
                i2 = p2[idx]
                if i1 != i2:
                    return self.fail("Sequence '{}': character {}: expected {} but found {}".format(label,
                        idx, i1, i2))
        return TestRunner.PASS

    def test_set_value_symmetric_difference(self):

        # v1 = [2, 2, 2, 3, 4, 8, 12, 14, 3, 99, 39, 12]
        # v2 = [2, 2, 3, 8, 12, 8, 8, 8, 3, 17, 18, 20, 21, 33]
        # c1 = collections.Counter(v1)
        # c2 = collections.Counter(v2)
        # result = (c1 - c2) + (c2 - c1)
        # unmatched = list(result.elements())
        # print len(unmatched)
        # print " ".join(["{}".format(x) for x in sorted(unmatched)])

        # data = sorted([self.rng.randint(0, 5000) for i in range(num_raw_points)])
        nreps = 5
        for rep in range(nreps):
            v1 = [self.rng.randint(1, 100) for i in range( self.rng.randint(100, 200) ) ]
            v2 = [self.rng.randint(1, 100) for i in range( self.rng.randint(100, 200) ) ]

            src_str = " ".join([str(i) for i in v1]) + "\n" + " ".join([str(i) for i in v2])
            self.execute_test('set_value_symmetric_difference', None, src_str)
            if self.test_retcode != 0:
                return TestRunner.ERROR

            c1 = collections.Counter(v1)
            c2 = collections.Counter(v2)
            counted_matched = c1 & c2
            exp_matched = sorted(list(counted_matched.elements()))
            counted_diffs = (c1 - c2) + (c2 - c1)
            exp_unmatched = sorted(list(counted_diffs.elements()))
            exp_diff = len(exp_unmatched)

            rows = self.test_stdout.split("\n")
            obs_diff = int(rows[0])
            obs_matched = [int(i) for i in rows[1].split()]
            obs_unmatched = [int(i) for i in rows[2].split()]
            obs_diff2 = int(rows[3])

            debug_str = []
            debug_str.append("  Expected matched ({}):\n  {}".format( len(exp_matched)   , exp_matched,   ))
            debug_str.append("  Observed matched ({}):\n  {}".format( len(obs_matched)   , obs_matched,   ))
            debug_str.append("Expected unmatched ({}):\n  {}".format( len(exp_unmatched) , exp_unmatched, ))
            debug_str.append("Observed unmatched ({}):\n  {}".format( len(obs_unmatched) , obs_unmatched, ))
            debug_str.append("Expecting reported difference of {}, but found {}".format(exp_diff, obs_diff))
            debug_str = "\n".join(debug_str)

            if obs_diff != exp_diff:
                return self.fail("Incorrect reported difference:\n{}".format(debug_str))
            if obs_matched != exp_matched:
                return self.fail("Incorrect observed matched:\n{}".format(debug_str))
            if obs_unmatched != exp_unmatched:
                return self.fail("Incorrect observed unmatched:\n{}".format(debug_str))
            if obs_diff != obs_diff2:
                return self.fail("Internal counting error: {} vs. {}".format(obs_diff, obs_diff2))

        return TestRunner.PASS

    def test_profile_generation(self):
        # a = []
        # b = []
        # for raw_data_size in range(2, 1000):
        #     raw_data = range(1, raw_data_size+1)
        #     rds = len(raw_data)
        #     self.execute_test('profile_generation', [1000], " ".join([str(d) for d in raw_data]))
        #     out = [x for x in self.test_stdout.split("\n") if x]
        #     diff = len(out) - 1000
        #     if diff:
        #         a.append( (raw_data_size, len(out), diff) )
        #     else:
        #         b.append( (raw_data_size, len(out), diff) )
        # a.sort(key = lambda x: abs(x[2]))
        # for i in a:
        #     print i
        # print "--"
        # print len(a), len(b)

        num_raw_points = 1000
        num_profile_points = num_raw_points * 10
        # data = sorted([self.rng.randint(0, 5000) for i in range(num_raw_points)])
        data = range(1, num_raw_points)
        self.execute_test('profile_generation', [num_profile_points], " ".join([str(d) for d in data]))
        if self.test_retcode != 0:
            return TestRunner.ERROR
        observed = [float(i) for i in self.test_stdout.split("\n") if i]
        if len(observed) != num_profile_points:
            return self.fail("Expecting {} elements in profile but found {}: {}".format(
                    num_profile_points,
                    len(observed),
                    observed))
        data_idx = 0
        groups = {}
        for ov in observed:
            d = data[data_idx]
            while (ov > d):
                data_idx += 1
                d = data[data_idx]
            try:
                groups[d].append(ov)
            except KeyError:
                groups[d] = [ov]

        # make sure original points found in profile
        s1 = set(data)
        s2 = set(groups.keys())
        if s1 != s2:
            return self.fail("Not all raw data points represented in profile: {} (difference = {})".format(
                    s2,
                    s1.difference(s2)))

        # ensure bin sizes are roughly equal
        raw_data_value_min = None
        num_interpolated_points_min = float("inf")
        raw_data_value_max = None
        num_interpolated_points_max = 0
        for key, group in groups.items():
            val = len(group)
            if val < num_interpolated_points_min:
                raw_data_value_min = key
                num_interpolated_points_min = val
            if val > num_interpolated_points_max:
                raw_data_value_max = key
                num_interpolated_points_max = val
        if abs(num_interpolated_points_min - num_interpolated_points_max) > 1:
            return self.fail("Difference between minimum and maximum number of points interpolated between raw data points > 1: Raw data value '{}' results in '{}' interpolated points, while raw data value '{}' results in '{}' interpolated points".format(
                raw_data_value_min, num_interpolated_points_min,
                raw_data_value_max, num_interpolated_points_max))


        # bin_sizes = [len(v) for v in groups.values()]
        # min_bin_size = min(bin_sizes)
        # max_bin_size = max(bin_sizes)
        # print groups
        # print max_bin_size, min_bin_size

        return TestRunner.PASS

    def check_profile_distance(self, num_profiles, num_subprofiles, profile_size, max_data_val=1e4):

        profiles = []

        for i in range(num_profiles):
            subprofiles = []
            for j in range(num_subprofiles):
                subprofile = [self.rng.uniform(0, max_data_val) for j in range(profile_size)]
                noise = self.rng.uniform(-2, 10)
                while (noise > 0):
                    subprofile.insert(0, self.rng.uniform(0, max_data_val))
                    noise -= 1
                subprofiles.append(subprofile)
            profiles.append(subprofiles)

        input_data = [ ]
        for profile in profiles:
            for subprofile in profile:
                s = " ".join(["{:40.20f}".format(d) for d in subprofile])
                input_data.append(s)

        expected = []
        for profile_idx, profile1 in enumerate(profiles[:-1]):
            for profile2 in profiles[profile_idx+1:]:
                ss = 0
                for subprofile_idx, subprofile1 in enumerate(profile1):
                    subprofile2 = profile2[subprofile_idx]
                    if len(subprofile1) > len(subprofile2):
                        data1 = subprofile1[len(subprofile1) - len(subprofile2):]
                        data2 = subprofile2
                    elif len(subprofile2) > len(subprofile1):
                        data1 = subprofile1
                        data2 = subprofile2[len(subprofile2) - len(subprofile1):]
                    else:
                        data1 = subprofile1
                        data2 = subprofile2

                    assert len(data1) == len(data2)
                    assert len(data1) >= profile_size

                    for eidx in range(len(data1)):
                        v = pow(data1[eidx] - data2[eidx], 2)
                        ss += v
                dist = math.sqrt(ss)
                expected.append(dist)

        self.execute_test('profile_distances', [num_subprofiles], "\n".join(input_data))
        if self.test_retcode != 0:
            # print "\n***********************"
            # print self.test_stdout
            # print "---"
            # print self.test_stderr
            # print "***********************\n"
            return TestRunner.ERROR
        observed = [float(i) for i in self.test_stdout.split("\n") if i]

        if len(expected) != len(observed):
            # print "\n***********************"
            # print len(input_data)
            # print "---"
            # print self.test_stdout
            # print "---"
            # print self.test_stderr
            # print "***********************\n"
            return self.fail("Expecting {} elements in results but found {}".format(
                len(expected), len(observed)))

        for idx, expv in enumerate(expected):
            obsv = observed[idx]
            if not self.is_almost_equal(obsv, expv, 1e-3):
                return self.fail("Distance result #{}: expecting {} but found {}".format(
                    idx+1,
                    expv,
                    obsv))
        return TestRunner.PASS

    def test_simple_profile_distance(self):
        return self.check_profile_distance(10, 1, 1000)

    def test_multisubprofile_profile_distance(self):
        return self.check_profile_distance(10, 10, 1000)

    def compare_tree_splits_distances(self, treefile_basename, schema="nexus"):
        treefile = os.path.join(self.data_dir, "general", treefile_basename)
        self.execute_test("tree_splits_distances",
                [treefile, schema])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        trees = dendropy.TreeList.get_from_path(treefile, schema)
        obs_rows = self.test_stdout.split("\n")
        ridx = 0
        for tidx1, tree1 in enumerate(trees[:-1]):
            for tidx2, tree2 in enumerate(trees[tidx1:]):
                fp, fn = tree1.false_positives_and_negatives(tree2)
                sd = fp + fn
                ed = tree1.euclidean_distance(tree2)
                urfd = tree1.robinson_foulds_distance(tree2)
                labels = ["false negatives", "false positives", "symmetric difference", "Euclidean distance", "unweighted Robinson-Foulds distance"]
                exp_cols = [fn, fp, sd, ed, urfd]
                obs_cols = [float(v) for v in obs_rows[ridx].split("\t")]
                for cidx, expv in enumerate(exp_cols):
                    obsv = obs_cols[cidx]
                    if not self.is_almost_equal(expv, obsv):
                        return self.fail("Comparison {} vs. {}: expecting {} of {} but found {} (difference = {}); <{}> vs. <{}>".format(
                        tidx1 + 1,
                        tidx2 + 2,
                        labels[cidx],
                        expv,
                        obsv,
                        expv - obsv,
                        ", ".join([str(v) for v in exp_cols]),
                        ", ".join([str(v) for v in obs_cols]),
                        ))
                ridx += 1
        return TestRunner.PASS

    def test_tree_splits_distance(self):
        # return self.compare_tree_splits_distances("/Users/jeet/Scratch/x.tre", "newick")
        return self.compare_tree_splits_distances("pythonidae.reference-trees.nexus", "nexus")

    def test_tree_subtree_leaf_set_sizes(self):
        treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
        self.execute_test("subtree_leaf_set_sizes",
                [treefile, "newick"])
        if self.test_retcode != 0:
            return TestRunner.ERROR
        tree = dendropy.Tree.get_from_path(treefile, "newick")
        internal_nodes = tree.internal_nodes()
        exp_subtree_leaf_set_sizes = {}
        for nd in internal_nodes:
            leaf_count = 0
            for leaf in nd.leaf_iter():
                leaf_count += 1
            exp_subtree_leaf_set_sizes[nd.label] = leaf_count
        obs_subtree_leaf_set_sizes = {}
        obs_rows = self.test_stdout.split("\n")
        for row in obs_rows:
            if not row:
                obs_rows.remove(row)
                continue
            obs_cols = row.split("\t")
            label = obs_cols[0]
            sz1 = int(obs_cols[1])
            sz2 = int(obs_cols[2])
            if label not in exp_subtree_leaf_set_sizes:
                return self.fail("Unexpected node: '{}'".format(label))
            exp_sz = exp_subtree_leaf_set_sizes[label]
            if sz1 != sz2:
                return self.fail("Internal check mismatch in leaf set sizes: '{}': {} vs. {}".format(label, sz1, sz2))
            if exp_sz != sz1:
                return self.fail("External check mismatch in leaf set sizes: '{}': expecting {} but found {}".format(label, exp_sz, sz1))
        if len(obs_rows) != len(exp_subtree_leaf_set_sizes):
            return self.fail("Expecting {} leaf set sizes but found {}".format(len(exp_subtree_leaf_set_sizes), len(obs_rows)))
        return TestRunner.PASS

    # def test_tree_subtree_clade_sizes(self):
    #     treefile = os.path.join(self.data_dir, "general", "pythonidae.postorder.newick")
    #     # treefile = os.path.join(self.data_dir, "general", os.path.expanduser("~/Scratch/x.tre"))
    #     self.execute_test("subtree_clade_sizes",
    #             [treefile, "newick"])
    #     if self.test_retcode != 0:
    #         return TestRunner.ERROR
    #     tree = dendropy.Tree.get_from_path(treefile, "newick")
    #     internal_nodes = tree.internal_nodes()
    #     exp_subtree_clade_sizes = {}
    #     for nd in internal_nodes:
    #         nd_count = -1 # adjust for top node being counted in iteration
    #         for nd in nd.postorder_iter():
    #             nd_count += 1
    #         exp_subtree_clade_sizes[nd.label] = nd_count
    #     obs_subtree_clade_sizes = {}
    #     obs_rows = self.test_stdout.split("\n")
    #     for row in obs_rows:
    #         if not row:
    #             obs_rows.remove(row)
    #             continue
    #         obs_cols = row.split("\t")
    #         label = obs_cols[0]
    #         sz1 = int(obs_cols[1])
    #         if label not in exp_subtree_clade_sizes:
    #             return self.fail("Unexpected node: '{}'".format(label))
    #         exp_sz = exp_subtree_clade_sizes[label]
    #         if exp_sz != sz1:
    #             return self.fail("External check mismatch in clade size: '{}': expecting {} but found {}".format(label, exp_sz, sz1))
    #             # print
    #             # print self.test_stdout
    #             # print tree.as_ascii_plot(show_internal_node_labels=True)
    #     if len(obs_rows) != len(exp_subtree_clade_sizes):
    #         return self.fail("Expecting {} clade sizes but found {}".format(len(exp_subtree_clade_sizes), len(obs_rows)))
    #     return TestRunner.PASS

    def run(self):
        tests_to_run = []
        isolated_tests_to_run = []
        for name, value in inspect.getmembers(self, callable):
            if name.startswith("test"):
                tests_to_run.append((name, value))
            elif name.startswith("isolate_test"):
                isolated_tests_to_run.append((name, value))
        if isolated_tests_to_run:
            tests_to_run = isolated_tests_to_run
        passes = []
        fails = []
        errors = []
        self.cout("Running with random seed: {}\n".format(self.random_seed))
        for test_idx, (test_name, test_call) in enumerate(tests_to_run):
            self.cout("@turquoise@<<{: 4d}/{:<4d}>>: {}: ".format(test_idx+1, len(tests_to_run), test_name))
            result = test_call()
            if result == TestRunner.PASS:
                self.cout("@green@<<PASS>>\n")
                passes.append(test_name)
                if self.test_pass_message and self.verbosity > 3:
                    self.cout("         :   - {}\n".format(self.test_pass_message))
            elif result == TestRunner.FAIL:
                self.cout("@fuchsia@<<FAIL>>\n")
                fails.append(test_name)
                if self.test_fail_message:
                    self.cout("         :   - {}\n".format(self.test_fail_message))
                if self.break_on_fail:
                    self.summarize(passes, fails, errors)
                    return
            elif result == TestRunner.ERROR:
                self.cout("@red@<<ERROR>>\n")
                self.cout("         @red@<<:        Executed:>> {}\n".format(self.test_command))
                self.cout("         @red@<<:     Return Code:>> {}\n".format(self.test_retcode))
                self.cout("         @red@<<: Standard Output:>> {}\n".format(self.test_stdout))
                self.cout("         @red@<<:  Standard Error:>> {}\n".format(self.test_stderr))
                errors.append(test_name)
                if self.break_on_fail:
                    self.summarize(passes, fails, errors)
                    return
            else:
                raise Exception("Test did not return a valid result")
        self.summarize(passes, fails, errors)

    def summarize(self, passes, fails, errors):
        self.cout("\n--\nTests completed.\n")
        self.cout("@turquoise@<<{}>> tests run with @green@<<{}>> successes, @fuchsia@<<{}>> failures and @red@<<{}>> errors.\n".format(
                len(passes)+len(fails)+len(errors),
                len(passes),
                len(fails),
                len(errors)))

def main():
    """
    Main CLI handler.
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbosity",
            default=1,
            help="control messaging level")

    parser.add_argument("-b", "--break-on-fail",
            action="store_true",
            default=False,
            help="terminate tests after first failure")

    random_seed = random.randint(0, sys.maxint)
    parser.add_argument("-z", "--random-seed",
            default=random_seed,
            help="random seed")

    args = parser.parse_args()

    test_runner = TestRunner(args)
    test_runner.run()

if __name__ == '__main__':
    main()


