#!/usr/bin/env python3

import argparse
import pysam
import itertools

def parse_args():
    parser = argparse.ArgumentParser(description="Process BAM/SAM file and generate pairs based on mode.")
    parser.add_argument("-i", "--input", required=True, help="Input BAM/SAM file")
    parser.add_argument("-m", "--mode", required=True, choices=["greedy", "recursive"], help="Mode: greedy or recursive")
    return parser.parse_args()

def process_greedy(lst, readID):
    for x in itertools.combinations(lst, 2):
        tmp0 = x[0].split('#')
        tmp1 = x[1].split('#')
        print(f'{readID}\t{tmp0[0]}\t{tmp0[1]}\t{tmp1[0]}\t{tmp1[1]}\t{tmp0[2]}\t{tmp1[2]}\tUU')

def process_recursive(lst, readID):
    for i in range(len(lst) - 1):
        tmp0 = lst[i].split('#')
        tmp1 = lst[i + 1].split('#')
        print(f'{readID}\t{tmp0[0]}\t{tmp0[1]}\t{tmp1[0]}\t{tmp1[1]}\t{tmp0[2]}\t{tmp1[2]}\tUU')

def bam2pairs(samfile, mode):
    print('## pairs format v1.0')
    # 从 BAM/SAM 文件的 header 中获取染色体长度
    for rid, length in zip(samfile.references, samfile.lengths):
        print(f'#chromsize: {rid} {length}')
    print('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type')

    lst = []
    readID = ''
    for read in samfile:
        if not read.is_unmapped:
            current_readID = read.query_name.split(':')[0]
            if current_readID != readID:
                if readID != '':
                    if mode == 'greedy':
                        process_greedy(lst, readID)
                    elif mode == 'recursive':
                        process_recursive(lst, readID)
                lst = []
            readID = current_readID
            rid = read.reference_name
            start = read.reference_start + 1
            strand = '-' if read.is_reverse else '+'
            info = f'{rid}#{start}#{strand}'
            lst.append(info)

    # 处理最后一个readID
    if mode == 'greedy':
        process_greedy(lst, readID)
    elif mode == 'recursive':
        process_recursive(lst, readID)

def main():
    args = parse_args()

    samfile = pysam.AlignmentFile(args.input, "rb") if args.input.endswith('.bam') else pysam.AlignmentFile(args.input, "r")
    bam2pairs(samfile, args.mode)
    samfile.close()

if __name__ == "__main__":
    main()
