#! /usr/bin/env python

import sys
import random
import collections

rng = random.Random()

v1 = [rng.randint(1, 100) for i in range( rng.randint(100, 200) ) ]
v2 = [rng.randint(1, 100) for i in range( rng.randint(100, 200) ) ]
c1 = collections.Counter(v1)
c2 = collections.Counter(v2)
counted_matched = c1 & c2
matched = sorted(list(counted_matched.elements()))
counted_diffs = (c1 - c2) + (c2 - c1)
unmatched = sorted(list(counted_diffs.elements()))
diff = len(unmatched)

sys.stdout.write("{}\n".format(" ".join([str(i) for i in v1])))
sys.stdout.write("{}\n".format(" ".join([str(i) for i in v2])))
sys.stdout.write("{}\n".format(diff))
