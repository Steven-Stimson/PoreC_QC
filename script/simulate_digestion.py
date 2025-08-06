#!/usr/bin/env python

import argparse
import re
import os
import sys
import numpy as np
import pandas as pd

RE_cutsite = {
    "mboi": ["^GATC"],
    "dpnii": ["^GATC"],
    "bglii": ["A^GATCT"],
    "hindiii": ["A^AGCTT"]}

def find_re_sites(filename, sequences, offset):
    with open(filename, 'r') as infile:
        chr_id = None
        big_str = ""
        indices = []
        all_indices = []
        contig_names = []
        for line in infile:
            if line.startswith(">"):
                if chr_id is not None:
                    for rs in range(len(sequences)):
                        pattern = f"(?={sequences[rs].lower()})"
                        indices += [m.start() + offset[rs] 
                                   for m in re.finditer(pattern, big_str)]
                    indices.sort()
                    all_indices.append(indices.copy())
                    indices.clear()
                big_str = ""
                chr_id = line.split()[0][1:]
                contig_names.append(chr_id)
            else:
                big_str += line.lower().strip()
        # Process last chromosome
        for rs in range(len(sequences)):
            pattern = f"(?={sequences[rs].lower()})"
            indices += [m.start() + offset[rs] 
                       for m in re.finditer(pattern, big_str)]
        indices.sort()
        all_indices.append(indices)
    return contig_names, all_indices

def find_chromsome_lengths(filename):
    chr_names, lengths = [], []
    length = 0
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(">"):
                chr_names.append(line[1:].strip())
                if length != 0:
                    lengths.append(length)
                length = 0
            else:
                length += len(line.strip())
        lengths.append(length)
    return chr_names, np.array(lengths)

def replaceN(cs):
    if 'N' not in cs:
        return [cs]
    n_positions = [i for i, c in enumerate(cs) if c == 'N']
    permutations = []
    
    def _replace(cs, pos=0):
        if pos >= len(n_positions):
            permutations.append(cs)
            return
        current_pos = n_positions[pos]
        for nuc in ['A','C','G','T']:
            new_cs = cs[:current_pos] + nuc + cs[current_pos+1:]
            _replace(new_cs, pos+1)
    _replace(cs)
    return permutations

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="Input FASTA file")
    parser.add_argument('-r', '--restriction_sites', dest='res_sites', 
                        nargs='+', required=True,
                        help="Restriction sites format (e.g. 'A^AGCTT')")
    parser.add_argument('-e', '--efficiency', type=int, nargs='+', default=[100],
                        help="Efficiency values in percentages")
    parser.add_argument('-s', '--seed', type=int, 
                        help="Random seed for repeatable results")
    parser.add_argument('-o', '--output', default="fragment_lengths.csv",
                        help="Output TSV file path")
    
    args = parser.parse_args()

    # 参数验证
    if not all(0 <= e <= 100 for e in args.efficiency):
        print("Efficiency must be between 0-100", file=sys.stderr)
        sys.exit(1)
        
    np.random.seed(args.seed)

    # 处理酶切位点
    cutsites = []
    for s in args.res_sites:
        cutsites.extend(s.split(','))
        
    sequences, offsets = [], []
    for cs in cutsites:
        if cs.lower() in RE_cutsite:
            cseq = RE_cutsite[cs.lower()][0]
        else:
            cseq = cs
            
        offpos = cseq.find('^')
        if offpos == -1:
            print(f"Missing '^' in {cs}", file=sys.stderr)
            exit(1)
            
        # 处理无效字符
        invalid = [c for c in cseq if c not in ['A','T','G','C','N','^']]
        if invalid:
            print(f"Invalid characters {invalid} in {cs}", file=sys.stderr)
            exit(1)
            
        sequences.append(cseq.replace('^',''))
        offsets.append(offpos)

    # 展开包含N的位点
    expanded_sequences = []
    expanded_offsets = []
    for seq, off in zip(sequences, offsets):
        expanded = replaceN(seq)
        expanded_sequences.extend(expanded)
        expanded_offsets.extend([off]*len(expanded))

    # 获取所有可能的切割位置
    contig_names, all_sites = find_re_sites(
        args.fastafile, expanded_sequences, expanded_offsets)
    chr_names, lengths = find_chromsome_lengths(args.fastafile)

    # 数据收集
    results = []
    for eff in args.efficiency:
        print(f"Processing {eff}% digestion...")
        np.random.seed(args.seed)  # 每次效率重新设置种子以保持相关性
        
        # 处理每个染色体
        for chrom_idx in range(len(contig_names)):
            chr_name = contig_names[chrom_idx]
            total_length = lengths[chrom_idx]
            sites = all_sites[chrom_idx].copy()  # 获取原始位点
            
            # 应用效率过滤
            if eff < 100 and sites:
                prob = eff / 100.0
                mask = np.random.rand(len(sites)) < prob
                filtered_sites = np.array(sites)[mask].tolist()
            else:
                filtered_sites = sites
                
            # 计算片段长度
            starts = [0] + filtered_sites
            ends = filtered_sites + [total_length]
            for s, e in zip(starts, ends[1:]):
                if e > s:
                    results.append({
                        'Type': f"cut{eff}%",
                        'Length': e - s
                    })

    # 创建DataFrame并输出
    df = pd.DataFrame(results)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Output saved to {args.output}")
