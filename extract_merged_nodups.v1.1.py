#!/usr/bin/env python3

import sys
import itertools
import random
import pickle
from collections import defaultdict
import pysam
import argparse 

parser = argparse.ArgumentParser(
                    prog='The process of generating intermediate files for resolution caltulate',
                    description='The process of generating intermediate files for resolution caltulate')

parser.add_argument("--bam", type=str, required=True, help="List of BAM files sep by ','")
parser.add_argument("--mapq_threshold", type=int, default=10, help="MAPQ threshold for filtering")
# 解析参数
args = parser.parse_args()
# 使用参数
filepath = args.bam.rstrip().split(',')
mapq_threshold = args.mapq_threshold


def process_bam_file(inFile,mapq_th=0):
    if inFile.endswith('.bam'):
        samfile = pysam.AlignmentFile(inFile, "rb")
    else:
        samfile = pysam.AlignmentFile(inFile, "r")

    readID = ''
    lst = []

    for read in samfile:
        if not read.is_unmapped and read.mapping_quality>=mapq_th:
            if read.query_name.split(':')[0] == readID:
                readID = read.query_name.split(':')[0]
                rid = read.reference_name
                start = read.reference_start + 1
                if read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'
                mapq = read.mapping_quality
                info = '%s#%s#%s#%s' % (rid, start, strand, mapq)
                lst.append(info)
            else:
                if readID != '':
                    for x in itertools.combinations(lst, 2):
                        tmp0 = x[0].split('#')
                        tmp1 = x[1].split('#')
                        if tmp0[2] == "+":
                            strand_out0 = 0
                        else:
                            strand_out0 = 16
                        if tmp1[2] == "+":
                            strand_out1 = 0
                        else:
                            strand_out1 = 16
                        print(f'{strand_out0}\t{tmp0[0]}\t{tmp0[1]}\t1\t{strand_out1}\t{tmp1[0]}\t{tmp1[1]}\t2\t{tmp0[3]}\t.\t.\t{tmp1[3]}\t.\t.\t.\t.\t')
                lst = []
                readID = read.query_name.split(':')[0]
                rid = read.reference_name
                start = read.reference_start + 1
                if read.is_reverse:
                    strand = '-'
                else:
                    strand = '+'
                mapq = read.mapping_quality
                info = '%s#%s#%s#%s' % (rid, start, strand, mapq)
                lst.append(info)
            
            readID = read.query_name.split(':')[0]

    for x in itertools.combinations(lst, 2):
        tmp0 = x[0].split('#')
        tmp1 = x[1].split('#')
        if tmp0[2] == "+":
            strand_out0 = 0
        else:
            strand_out0 = 16
        if tmp1[2] == "+":
            strand_out1 = 0
        else:
            strand_out1 = 16
        print(f'{strand_out0}\t{tmp0[0]}\t{tmp0[1]}\t1\t{strand_out1}\t{tmp1[0]}\t{tmp1[1]}\t2\t{tmp0[3]}\t.\t.\t{tmp1[3]}\t.\t.\t.\t.\t')


for bam_file in filepath:
    process_bam_file(bam_file,mapq_threshold)
