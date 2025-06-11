#!/usr/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--order_type_file', required=True, type=str, help='QC.py output <04.xx.order_type.split.txt>')
parser.add_argument('--read_type_result_file', required=True, type=str, help='QC.py output <01.xx.read_type.txt>')
parser.add_argument("--out_prefix", type=str, required=True, help="Prefix for the output files")
args = parser.parse_args()
filename = args.order_type_file
read_type = args.read_type_result_file
out_prefix = args.out_prefix

df = pd.read_csv(filename, sep='\t', header=0)
df_read_type = pd.read_csv(read_type, sep='\t', header=0)

df2 = df.groupby("order").agg({'intra_read_num': 'sum', 'inter_read_num': 'sum'}).apply(lambda x: round(x *100 / sum(df_read_type["filter_mapping_read"]),2))

outfile = '06.{}.order_type.merge.txt'.format(out_prefix)
df2.to_csv(outfile, sep="\t")
