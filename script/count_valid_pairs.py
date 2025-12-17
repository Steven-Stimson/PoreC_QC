#!/usr/bin/env python3

import sys
import pandas as pd
import pyranges as pr
import gzip

if len(sys.argv) != 2:
    sys.exit('python3 %s <add.pairs.gz>' % (sys.argv[0]))
    
inFile = sys.argv[1]

total = 0
valid = 0
intra = 0
inter = 0

def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

with open_maybe_gz(inFile) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith('#'):
            continue
        total += 1
        tmp = line.split('\t')
        chrA = tmp[1]
        chrB = tmp[3]
        regionA = tmp[8]
        regionB = tmp[9]
        if regionA != regionB:
            valid += 1
            if chrA == chrB:
                intra += 1
            else:
                inter += 1
                
print('total\t%i' % (total))
print('valid\t%i' % (valid))
print('intra\t%i' % (intra))
print('inter\t%i' % (inter))
