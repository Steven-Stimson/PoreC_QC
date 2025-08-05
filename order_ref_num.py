#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
        '--input',
        '-i',
        type= str,
        help='order_read_type.stat'
)

parser.add_argument(
        '--ID',
        '-o',
        type= str,
        help='prefix out output'
)

args = parser.parse_args()


df = pd.read_csv(args.input, sep='\t')


# 设置order的排序顺序
order_categories = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '>=10']
df['order'] = pd.Categorical(df['order'], categories=order_categories, ordered=True)

# 处理ref_num：将>=5转为数值5
df['ref_num'] = df['ref_num'].replace('>=5', 5).astype(int)

# 创建透视表
pivot_df = df.pivot(index='order', columns='ref_num', values='ratio').fillna(0)

# 排列ref_num顺序并生成颜色
cols = sorted(pivot_df.columns)
pivot_df = pivot_df[cols]
colors = plt.cm.YlGn(np.linspace(1, 0.2, len(cols)))

# 绘制图表
fig, ax = plt.subplots(figsize=(12, 6))
pivot_df.plot.bar(stacked=True, ax=ax, color=colors, width=0.8, edgecolor='black')

# 设置标签和图例
ax.set_xlabel('Order')
ax.set_ylabel('Ratio')
ax.set_title(args.ID)
handles, labels = ax.get_legend_handles_labels()
ax.tick_params(axis='x', rotation=0)


# 处理图例标签：将5显示为'>=5'
new_labels = [str(l) if int(l)<5 else '>=5' for l in labels]
ax.legend(handles, new_labels, title='Ref_num', 
		bbox_to_anchor=(1.05, 1), loc='upper left')

# 调整布局并保存
plt.tight_layout()
plt.savefig(f'{args.ID}.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.savefig(f'{args.ID}.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig(f'{args.ID}.svg', format='svg', bbox_inches='tight')
