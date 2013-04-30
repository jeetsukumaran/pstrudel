#! /usr/bin/env python

import sys
import dendropy
from dendropy import treecalc

trees = dendropy.TreeList.get_from_path(sys.argv[1], sys.argv[2])
for tree in trees:
    pdm = treecalc.PatristicDistanceMatrix(tree)
    for tidx1 in range(len(tree.taxon_set)-1):
        for tidx2 in range(tidx1+1, len(tree.taxon_set)):
            taxon1 = tree.taxon_set[tidx1]
            taxon2 = tree.taxon_set[tidx2]
            exp_steps = pdm.path_edge_count(taxon1, taxon2)
            exp_weight = pdm(taxon1, taxon2)
            sys.stdout.write("{}\t{}\t{}\t{:0.40f}\n".format(taxon1.label, taxon2.label, exp_steps, exp_weight))
