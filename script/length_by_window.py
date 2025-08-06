#!/usr/bin/env python
import argparse
import pandas as pd
import numpy as np

def create_bins(lengths, window_size):
    """Create bins for length distribution"""
    min_val = 0
    max_val = np.ceil(lengths.max() / window_size) * window_size
    bins = np.arange(min_val, max_val+1, window_size)
    return bins

def bin_length(length, bins):
    """Get the bin label for a length value"""
    bin_label = f"{bins[np.digitize(length, bins)-1]}-{bins[np.digitize(length, bins)]}"
    return bin_label

def calculate_frequencies(df, window_size):
    """Calculate frequency statistics per type and bin"""
    all_results = []
    
    for type_group, group in df.groupby('Type'):
        bins = create_bins(group['Length'], window_size)
        group['Bin'] = pd.cut(
            group['Length'], 
            bins=bins, 
            right=False, 
            include_lowest=True
        ).apply(lambda x: f"{x.left:.0f}-{x.right:.0f}")

        counts = group['Bin'].value_counts().sort_index()
        total = len(group)
        
        for bin_range, count in counts.items():
            frequency = (count / total) * 100
            all_results.append({
                'Type': type_group,
                'Window': bin_range,
                'Frequency(%)': round(frequency, 2)
            })
    
    return pd.DataFrame(all_results)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Input TSV file path')
    parser.add_argument('-o', '--output', default='frequency_stats.tsv',
                        help='Output TSV file path')
    parser.add_argument('-w', '--window_size', type=int, required=True,
                        help='Window size for binning')
    
    args = parser.parse_args()
    
    # 读取TSV文件
    df = pd.read_csv(args.input_file, sep='\t')
    
    # 确保Length列是数值类型
    df['Length'] = df['Length'].astype(float)
    
    # 计算频率统计
    result_df = calculate_frequencies(df, args.window_size)
    
    # 保存结果
    result_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"Output saved to {args.output}")
    
if __name__ == "__main__":
    main()
