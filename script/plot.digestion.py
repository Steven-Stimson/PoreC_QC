#!/usr/bin/env python
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def filter_zero_windows(df):
    """过滤所有类型频率均为0的窗口"""
    filtered = df.groupby('Window_End').filter(
        lambda x: (x['Frequency(%)'] != 0).any()
    )
    return filtered

def plot_frequency(df, output, x_celling):
    """绘制分类型频率分布图（支持x轴范围控制）"""
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df,
        x='Window_End',
        y='Frequency(%)',
        hue='Type',
        marker='o',
        linewidth=2.5
    )
    plt.title('Fragment Length Frequency Distribution by Digestion Efficiency')
    plt.xlabel('Window(bp)')
    plt.ylabel('Frequency(%)')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(rotation=45)
    
    # 设置x轴范围（从0到指定最大值）
    plt.xlim(0, x_celling)
    
    plt.tight_layout()
    plt.savefig(f'{output}.png', format='png', dpi=300)
    plt.savefig(f'{output}.pdf', format='pdf', dpi=300)
    plt.savefig(f'{output}.svg', format='svg')

def main():
    parser = argparse.ArgumentParser(description="Fragment Frequency Visualization")
    parser.add_argument('input_file', help='输入统计TSV文件路径（包含Type/Window/Frequency列）')
    parser.add_argument('-o', '--output', default='digestion', 
                        help='prefix of output')
    parser.add_argument('--x_celling', type=int, default=20000, 
                        help='设置x轴的最大显示范围（默认20000）')
    
    args = parser.parse_args()
    
    df = pd.read_csv(args.input_file, sep='\t')
    df['Window_End'] = (
        df['Window']
        .str.extract(r'(\d+)-(\d+)')
        .iloc[:, 1]
        .astype(int)
    )
    
    filtered_df = filter_zero_windows(df)
    
    plot_frequency(filtered_df, args.output, args.x_celling)

if __name__ == "__main__":
    main()
