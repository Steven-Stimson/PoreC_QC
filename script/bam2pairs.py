#!/usr/bin/env python3

import sys
import itertools
import random
import pickle
from collections import defaultdict
import pysam

if len(sys.argv) != 2:
    sys.exit('python3 %s <inFile.bam>' % (sys.argv[0]))

inFile = sys.argv[1]

print('## pairs format v1.0')
dic = dict()
samfile = pysam.AlignmentFile(inFile, "rb") if inFile.endswith('.bam') else pysam.AlignmentFile(inFile, "r")

# 从 BAM/SAM 文件的 header 中获取染色体长度
for rid, length in zip(samfile.references, samfile.lengths):
    print(f'#chromsize: {rid} {length}')
print('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type')

lst = []

readID = ''
last_rid = ''
last_start = 0
last_end = 0
last_strand = ''
last_mapq = 0

for read in samfile:
    if not read.is_unmapped:
        if read.query_name.split(':')[0] == readID:
            readID = read.query_name.split(':')[0]
            rid = read.reference_name
            start = read.reference_start + 1
            end = read.reference_end
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            mapq = read.mapping_quality
            info = '%s#%s#%s' % (rid, start, strand)
            lst.append(info)
        else:
            if readID != '':
                for x in itertools.combinations(lst, 2):
                    tmp0 = x[0].split('#')
                    tmp1 = x[1].split('#')
                    print(f'{readID}\t{tmp0[0]}\t{tmp0[1]}\t{tmp1[0]}\t{tmp1[1]}\t{tmp0[2]}\t{tmp1[2]}\tUU')
            lst = []
            readID = read.query_name.split(':')[0]
            rid = read.reference_name
            start = read.reference_start
            end = read.reference_end
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            mapq = read.mapping_quality
            info = '%s#%s#%s' % (rid, start, strand)
            lst.append(info)

        readID = read.query_name.split(':')[0]

for x in itertools.combinations(lst, 2):
    tmp0 = x[0].split('#')
    tmp1 = x[1].split('#')
    print(f'{readID}\t{tmp0[0]}\t{tmp0[1]}\t{tmp1[0]}\t{tmp1[1]}\t{tmp0[2]}\t{tmp1[2]}\tUU')
