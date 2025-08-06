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

plt.figure(figsize=(12, 8))

types = df.columns[1:]
orders = df['order']

colors = plt.cm.YlGn(np.linspace(1, 0.2, len(orders)))


for j, type in enumerate(types):
	bottom = 0
	for i, order in enumerate(orders):
		plt.barh(j, df[type][i], left=bottom, height=0.5, color=colors[i], label=str(order) if j == 0 else "", edgecolor='black')
		bottom += df[type][i]

handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(
	by_label.values(), 
	by_label.keys(), 
	title='Reads Order', 
	loc='upper right',
	bbox_to_anchor=(1.2, 1),
	borderaxespad=0.
)

plt.xlabel('Percentage (%)')
plt.ylabel('Type')
plt.title(args.ID)

plt.yticks(range(len(types)), types)


for j, type in enumerate(types):
	total_percentage = df[type].sum()
	plt.text(0, j - 0.1, f'({total_percentage:.2f}%)', ha='right', va='center')

plt.tight_layout()
plt.savefig(f'{args.ID}.pdf',dpi=300, format='pdf', bbox_inches='tight')
plt.savefig(f'{args.ID}.png',dpi=300, format='png', bbox_inches='tight')
plt.savefig(f'{args.ID}.svg', format='svg', bbox_inches='tight')
