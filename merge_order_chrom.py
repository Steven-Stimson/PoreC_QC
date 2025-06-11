#!/usr/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--order_chrom_file', type=str, required=True, help='QC.py output <05.xx.order_chrom.split.txt>')
parser.add_argument("--out_prefix", type=str, required=True, help="Prefix for the output files")
args = parser.parse_args()
filename = args.order_chrom_file
out_prefix = args.out_prefix

df = pd.read_csv(filename, sep='\t', header=0)

#merge order_chrom_results
df1 = df.groupby(["order","chrom_num"]).agg({'read_num': 'sum'}).reset_index()
df2 = df1.groupby("order").agg({'read_num': 'sum'})
df3 = pd.merge( df1,df2, on="order", how="outer")
df3['read_ratio'] = round(df3['read_num_x'] / df3['read_num_y']*100,2)

outfile = '07.{}.order_chrom.merge.txt'.format(out_prefix)
df3.to_csv(outfile, index=False,sep="\t")
