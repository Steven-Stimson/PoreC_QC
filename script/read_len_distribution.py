#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description="绘制Read Length分布的计数图")
    parser.add_argument("-i", "--input", required=True, help="输入的TSV文件路径")
    parser.add_argument("--ID", required=True, help="输出文件的前缀")
    parser.add_argument("--window", type=int, default=500, help="窗口大小，默认为500")
    args = parser.parse_args()
    return args

def plot_read_len_count(input_tsv, out_prefix, window_size):
    """根据TSV文件绘制分类型计数分布图"""
    # 读取数据
    df = pd.read_csv(input_tsv, sep="\t")

    # 按窗口大小统计每个窗口内的read数量
    df['window'] = (df['read_len'] // window_size) * window_size
    df_grouped = df.groupby(['window', 'read_type']).size().reset_index(name='count')

    # 绘制计数分布图
    plt.figure(figsize=(12, 6))
    sns.lineplot(
        data=df_grouped,
        x='window',
        y='count',
        hue='read_type',
        marker='o',
        linewidth=2.5
    )
    plt.title('Read Length Distribution by Type')
    plt.xlabel('Window(bp)')
    plt.ylabel('Count')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig(f'{out_prefix}.png', format='png', dpi=300)
    plt.savefig(f'{out_prefix}.pdf', format='pdf', dpi=300)
    plt.savefig(f'{out_prefix}.svg', format='svg')
    plt.close()

def main():
    """主函数，解析参数并执行绘图"""
    args = parse_arguments()
    plot_read_len_count(args.input, args.ID, args.window)

if __name__ == "__main__":
    main()
