#! /usr/bin/env python

import sys
import dendropy

trees = dendropy.TreeList.get_from_path(sys.argv[1], sys.argv[2])
for tidx1 in range(len(trees)-1):
    for tidx2 in range(tidx1+1, len(trees)-1):
        fp, fn = trees[tidx1].false_positives_and_negatives(trees[tidx2])
        eed = trees[tidx1].euclidean_distance(trees[tidx2])
        sys.stdout.write("{}\t{}\t{}\t{}\n".format(fp, fn, fp + fn, eed))
