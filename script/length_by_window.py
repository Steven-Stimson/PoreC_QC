#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

def create_bins(lengths, window_size):
    min_val = 0
    max_val = np.ceil(lengths.max() / window_size) * window_size
    bins = np.arange(min_val, max_val+1, window_size)
    return bins

def calculate_frequencies(df, window_size):
    all_results = []
    
    # 确保Type按首次出现的顺序排列
    df['Type'] = pd.Categorical(df['Type'], categories=df['Type'].unique(), ordered=True)
    
    for type_group, group in df.groupby('Type', group_keys=False):
        bins = create_bins(group['Length'], window_size)
        group['Bin'] = pd.cut(
            group['Length'], 
            bins=bins, 
            right=False, 
            include_lowest=True
        ).apply(lambda x: f"{x.left:.0f}-{x.right:.0f}")
        
        # 获取Bin的首次出现顺序
        unique_bin_order = group['Bin'].unique()
        counts = group['Bin'].value_counts().reindex(unique_bin_order)
        
        total = len(group)
        
        for bin_range, count in counts.items():
            frequency = (count / total) * 100
            all_results.append({
                'Type': type_group,
                'Window': bin_range,
                'Frequency(%)': round(frequency, 2)
            })
    
    # 创建结果DataFrame并保持Type的顺序
    result_df = pd.DataFrame(all_results)
    # 确保Type顺序正确，并按Type和Bin的顺序排列
    result_df['Type'] = pd.Categorical(
        result_df['Type'], 
        categories=df['Type'].unique(), 
        ordered=True
    )
    result_df = result_df.sort_values(['Type', 'Window'])
    
    return result_df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Input TSV file path')
    parser.add_argument('-o', '--output', default='frequency_stats.tsv',
                        help='Output TSV file path')
    parser.add_argument('-w', '--window_size', type=int, required=True,
                        help='Window size for binning')
    
    args = parser.parse_args()
    
    df = pd.read_csv(args.input_file, sep='\t')
    df['Length'] = df['Length'].astype(float)
    
    result_df = calculate_frequencies(df, args.window_size)
    
    # 保存结果前，重置索引并保留原始顺序
    result_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"Output saved to {args.output}")

if __name__ == "__main__":
    main()
