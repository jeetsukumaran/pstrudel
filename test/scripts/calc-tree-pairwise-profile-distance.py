#! /usr/bin/env python

import sys
import math
import dendropy
from dendropy import treecalc

def euclidean_distance(v1, v2):
    ss = 0.0
    assert len(v1) == len(v2)
    for idx1 in range(len(v1)):
        ss += pow((v1[idx1] - v2[idx1]), 2)
    return math.sqrt(ss)

def main():
    trees = dendropy.TreeList.get_from_path(sys.argv[1], sys.argv[2])
    for tree in trees:
        pdm = treecalc.PatristicDistanceMatrix(tree)
        tree.wt_profile = sorted(pdm.distances())
        tree.uwt_profile = []
        for txidx1 in range(len(tree.taxon_set)-1):
            for txidx2 in range(txidx1+1, len(tree.taxon_set)):
                tree.uwt_profile.append(pdm.path_edge_count(tree.taxon_set[txidx1], tree.taxon_set[txidx2]))
        tree.uwt_profile = sorted(tree.uwt_profile)

    for tidx1 in range(len(trees)):
        for tidx2 in range(len(trees)):
            sys.stdout.write("{}\t{}\t{:.22f}\t{:.22f}\n".format(
                tidx1,
                tidx2,
                euclidean_distance(trees[tidx1].uwt_profile, trees[tidx2].uwt_profile),
                euclidean_distance(trees[tidx1].wt_profile, trees[tidx2].wt_profile)))

if __name__ == "__main__":
    main()



