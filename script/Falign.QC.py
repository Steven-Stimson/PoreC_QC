#!/usr/bin/env python3
import sys
import os
import pysam
from _collections import defaultdict
import argparse
import gc
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
        parser = argparse.ArgumentParser()
        parser.add_argument(
                '--frag_bam',
                '-i',
                type=str,
                help='paf file output by Falign'
        )

        parser.add_argument(
                '--Qthreshold',
                '-q',
                type=int,
                default=0,
                help='threshold of mapping quality,default 0',
        )
        parser.add_argument(
                '--ID',
                '-o',
                type=str,
                default='EntityID',
                help='prefix of output,default "EntityID"',
        )
    
        parser.add_argument(
                '--seqkit_stat',
                '-s',
                type=str,
                help='paf file output by Falign'
        )


        return parser.parse_args()


def calculate_n50(contig_lengths):
        # 检查输入列表是否为空
        if not contig_lengths:
                return 0  # 或者根据需求抛出异常

        # 1. 将contig长度降序排列
        sorted_contigs = sorted(contig_lengths, reverse=True)
    
        # 2. 计算总长度
        total_length = sum(sorted_contigs)
        half_total = total_length / 2.0
    
        # 3. 累计计算直到超过总长度的一半
        cumulative = 0
        for contig in sorted_contigs:
                cumulative += contig
                if cumulative >= half_total:
                        return contig  # 返回当前contig的长度作为N50
    
        return 0  # 如果在极端情况下所有contig加起来都不足总长的一半（几乎不可能）

def raw_data_info(seqkit_stat):
        stat = pd.read_table(seqkit_stat,sep='\s+')
        stat["num_seqs"] = stat["num_seqs"].str.replace(",", "").astype(int)
        stat["sum_len"] = stat["sum_len"].str.replace(",", "").astype(int)
        total_raw_readnum = sum(stat["num_seqs"])
        total_raw_baseslen = sum(stat["sum_len"])
        return total_raw_readnum, total_raw_baseslen

def MappingStat(bam, Qthreshold):


        samfile = pysam.AlignmentFile(bam, 'rb') if bam.endswith('.bam') else pysam.AlignmentFile(bam, "r")

        ReadsMap = defaultdict(lambda: {'sub_read_num': 0,
                                                                        'order': 0,
                                                                        'ref': defaultdict(int),
                                                                        'read_len': 0,
                                                                        'map_stat': list(),
                                                                        'read_type': str(),
                                                                        'ref_num': int(),
                                                                        'align': defaultdict(list)
                                                                        }
                                                                )
        references = samfile.references

        for read in samfile:
                main_read_id = (read.query_name).split(':')[0].rsplit('_', 1)[0]
                
                ReadsMap[main_read_id]['sub_read_num'] += 1
                ReadsMap[main_read_id]['read_len']=int(read.get_tag("ql"))
                
                if read.mapq > Qthreshold :
                        ReadsMap[main_read_id]['ref'][references[read.reference_id]] += 1
                        ReadsMap[main_read_id]['order'] += 1
                        ReadsMap[main_read_id]['align'][references[read.reference_id]].append([int(read.reference_start), int(read.reference_end)])
                else:
                        continue
                                                    
        for main_read_id in ReadsMap.keys():

                read_type = str()
                ref_num = int()

                if ReadsMap[main_read_id]['ref'].keys():
                        ref_num = len(ReadsMap[main_read_id]['ref'].keys())
                        max_map_ref = max(ReadsMap[main_read_id]['ref'].values())
                        total_map_ref = sum(ReadsMap[main_read_id]['ref'].values())
                        if ref_num == 1 and total_map_ref ==1 :
                                read_type = 'align-one'
                        elif max_map_ref > (total_map_ref - max_map_ref): #如果比上最多的那个参考比其余的加起来还要多
                                read_type = 'intra'
                        else:
                                read_type = 'inter'
                #else:
                #       read_type = 'unmap'
                #       ref_num = 0                                            
                ReadsMap[main_read_id]['read_type'] = read_type
                ReadsMap[main_read_id]['ref_num'] = ref_num                                            
                #ReadsMap[main_read_id]['read_n50'] = calculate_n50(ReadsMap[main_read_id]['read_len'])



        ReadsMapDF = pd.DataFrame.from_dict(ReadsMap, orient='index')

        del ReadsMap
        gc.collect()

        return ReadsMapDF

def output_read_type_stats(ReadsMapDF, total_raw_readnum, total_raw_baseslen, out_prefix):
        """输出四种read类型的统计信息"""
        # 计算read数目
        read_counts = ReadsMapDF['read_type'].value_counts()
        # 添加unmapped reads
        mapped_reads = len(ReadsMapDF)
        unmapped_reads = total_raw_readnum - mapped_reads
        read_counts['total'] = total_raw_readnum
        read_counts['unmapped'] = unmapped_reads
    
        # 计算base总数
        base_totals = ReadsMapDF.groupby('read_type')['read_len'].sum()
        mapped_bases = base_totals.sum()  # 已比对的总碱基数
        base_totals['total'] = total_raw_baseslen
        base_totals['unmapped'] = total_raw_baseslen - mapped_bases  # 直接计算未比对碱基
    
        # 转换为DataFrame并保存
        stats_df = pd.DataFrame({
                'read_count': read_counts,
                'base_total': base_totals
        }).fillna(0).astype(int)  # 填充NA并转为整数
    
        # 确保所有类型都存在
        for typ in ['total','align-one', 'intra', 'inter', 'unmapped']:
                if typ not in stats_df.index:
                        stats_df.loc[typ] = [0, 0]
    
        # 按固定顺序输出
        stats_df = stats_df.loc[['total','unmapped', 'align-one', 'intra', 'inter']]
    
        # 保存为TSV
        stats_df.to_csv(f"{out_prefix}.map_stats.tsv", sep="\t", float_format="%.0f")
    
        return stats_df



def order_read_type_read(ReadsMapDF, total_raw_readnum,outfile):
        df = ReadsMapDF.copy()
        #type_stat = df['read_type'].value_counts()

        def group_order(x):
                if x < 10:
                        return x
                else:
                        return '>=10'

        df['order'] = df['order'].apply(group_order)


        # 1. 计算每种order不同read_type的情况
        order_read_type_read = pd.crosstab(df['order'], df['read_type'])
        # 2. 添加unmapped reads统计
        mapped_reads = len(df)  # 已比对reads数
        unmapped_reads = total_raw_readnum - mapped_reads  # 未比对reads数
    
        if 0 not in order_read_type_read.index:
                order_read_type_read.loc[0] = [0, 0, 0]

    # 新增 'unmapped' 列，并给 order=0 的行赋值为 unmapped_reads
        order_read_type_read['unmapped'] = 0  # 先初始化为 0
        order_read_type_read.at[0, 'unmapped'] = unmapped_reads  # 修改 order=0 行的 'unmap' 值
    
        # 2. 除以reads总数得到不同read_type的比例
        order_read_type_read = order_read_type_read.divide(total_raw_readnum, axis=0) * 100
        # 3. 保留两位小数
        order_read_type_read = order_read_type_read.round(2)
        # 4. 调整顺序
        order_read_type_read = order_read_type_read[['unmapped', 'align-one', 'intra', 'inter']]
        # 5. 输出为tsv格式
        order_read_type_read.to_csv(outfile, sep='\t', index=True)

        return order_read_type_read

def order_read_type_base(ReadsMapDF, total_raw_baseslen,outfile):
        df = ReadsMapDF.copy()
        #type_stat = df['read_type'].value_counts()

        def group_order(x):
                if x < 10:
                        return x
                else:
                        return '>=10'

        df['order'] = df['order'].apply(group_order)


        # 1. 计算每种order不同read_type的情况
        order_read_type_base = df.pivot_table(
                index='order',        # 行索引
                columns='read_type',  # 列索引
                values='read_len',    # 要聚合的值
                aggfunc='sum',        # 聚合方式：求read_len和
                fill_value=0          # 缺失值填充为 0
        )

        # 2. 添加unmapped reads统计
        mapped_base = sum(df["read_len"])  # 已比对base数
        unmapped_base = total_raw_baseslen - mapped_base  # 未比对base数
    
        if 0 not in order_read_type_base.index:
                order_read_type_base.loc[0] = [0, 0, 0]

    # 新增 'unmapped' 列，并给 order=0 的行赋值为 unmapped_reads
        order_read_type_base['unmapped'] = 0  # 先初始化为 0
        order_read_type_base.at[0, 'unmapped'] = unmapped_base  # 修改 order=0 行的 'unmap' 值
    
        # 2. 除以reads总数得到不同read_type的比例
        order_read_type_base = order_read_type_base.divide(total_raw_baseslen, axis=0) * 100
        # 3. 保留两位小数
        order_read_type_base = order_read_type_base.round(2)
        # 4. 调整顺序
        order_read_type_base = order_read_type_base[['unmapped', 'align-one', 'intra', 'inter']]
        # 5. 输出为tsv格式
        order_read_type_base.to_csv(outfile, sep='\t', index=True)

        return order_read_type_base



def order_ref_num_read(ReadsMapDF, outfile):
        df = ReadsMapDF.copy()

        def group_order(x):
                if x < 10:
                        return x
                else:
                        return '>=10'

        def group_ref_num(x):
                if x < 5:
                        return x
                else:
                        return '>=5'

        df['order'] = df['order'].apply(group_order)
        df['ref_num'] = df['ref_num'].apply(group_ref_num)

        # 1. 统计每个order-ref组合的出现次数（count列）
        order_ref_num_read = df.groupby(['order', 'ref_num']).size().reset_index(name='count')
        # 2. 计算每个order的总数（total列）
        order_ref_num_read['total'] = order_ref_num_read.groupby('order')['count'].transform('sum')
        # 3. 计算占比（ratio列）并保留四位小数
        order_ref_num_read['ratio'] = (order_ref_num_read['count'] / order_ref_num_read['total']).round(2)
        # 4. 调整列顺序为：order → ref_num → count → total → ratio
        order_ref_num_read = order_ref_num_read[['order', 'ref_num', 'count', 'total', 'ratio']]
        # 5. 输出为tsv格式
        order_ref_num_read.to_csv(outfile, sep='\t', index=False)


        return order_ref_num_read


def order_ref_num_base(ReadsMapDF, outfile):
        df = ReadsMapDF.copy()

        def group_order(x):
                if x < 10:
                        return x
                else:
                        return '>=10'

        def group_ref_num(x):
                if x < 5:
                        return x
                else:
                        return '>=5'

        df['order'] = df['order'].apply(group_order)
        df['ref_num'] = df['ref_num'].apply(group_ref_num)

        # 1. 统计每个order-ref组合的出现次数（count列）
        order_ref_num_base = df.groupby(['order', 'ref_num'])['read_len'].sum().reset_index(name='sum_read_len')
        # 2. 计算每个order的总数（total列）
        order_ref_num_base['total'] = order_ref_num_base.groupby('order')['sum_read_len'].transform('sum')
        # 3. 计算占比（ratio列）并保留四位小数
        order_ref_num_base['ratio'] = (order_ref_num_base['sum_read_len'] / order_ref_num_base['total']).round(2)
        # 4. 调整列顺序为：order → ref_num → count → total → ratio
        order_ref_num_base = order_ref_num_base[['order', 'ref_num', 'sum_read_len', 'total', 'ratio']]
        # 5. 输出为tsv格式
        order_ref_num_base.to_csv(outfile, sep='\t', index=False)

        return order_ref_num_base


def read_len_order(ReadsMapDF, outfile):
        pass



def map_length(ReadsMapDF, outfile):
        df = ReadsMapDF.copy()
        df = df[['map_stat', 'read_len']]
        df = df.explode(['map_stat', 'read_len']).reset_index(drop=True)
        df.to_csv(outfile, sep='\t', index=False)



def cutter_gap(ReadsMapDF, outfile):
        df = ReadsMapDF.copy()
        cutter_gap = []
        for _, row in df.iterrows():
                #order = row['order']
                #if order == 0: continue
                align_dict = row['align']

                for chrom, intervals in align_dict.items():
                        # 按比对起点排序
                        sorted_intervals = sorted(intervals, key=lambda x: x[0])
        
                # 计算相邻区间的gap
                for i in range(len(sorted_intervals)-1):
                        prev_end = sorted_intervals[i][1]
                        next_start = sorted_intervals[i+1][0]
                        gap = abs(next_start - prev_end)
            
                        cutter_gap.append({
                                'gap': gap,
                                #'order': order
                        })

        # 创建结果DataFrame
        cutter_gap_df = pd.DataFrame(cutter_gap)

        cutter_gap_df.to_csv(outfile, sep='\t', index=False)

def plot_read_len_distribution(ReadsMapDF, out_prefix):
    """绘制三类read的长度分布图"""
    df = ReadsMapDF.copy()
    
    if df.empty:
        print("Warning: No reads for plotting distribution")
        return

    # 绘制直方图
    plt.figure(figsize=(12, 6))
    for typ in ['align-one', 'intra', 'inter']:
        subset = df[df['read_type'] == typ]
        if not subset.empty:
            plt.hist(subset['read_len'], bins=50, alpha=0.5, label=typ)
    plt.xlabel('Read Length')
    plt.ylabel('Count')
    plt.xlim(0, 20000)
    plt.title('Read Length Distribution')
    plt.legend()
    plt.tight_layout()
    
    # 保存直方图
    plt.savefig(f"{out_prefix}.read_len_distribution.histogram.png", format='png', dpi=300)
    plt.savefig(f"{out_prefix}.read_len_distribution.histogram.pdf", format='pdf', dpi=300)
    plt.savefig(f"{out_prefix}.read_len_distribution.histogram.svg", format='svg')

    plt.close()

    # 绘制箱线图
    plt.figure(figsize=(12, 6))
    data_to_plot = [df[df['read_type'] == typ]['read_len'] for typ in ['align-one', 'intra', 'inter']]
    plt.boxplot(data_to_plot, labels=['align-one', 'intra', 'inter'])
    plt.xlabel('Read Type')
    plt.ylabel('Read Length')
    plt.ylim(0, 20000) 
    plt.title('Read Length by Type')
    plt.tight_layout() 
    # 保存箱线图
    plt.savefig(f"{out_prefix}.read_len_distribution.boxplot.png", format='png', dpi=300)
    plt.savefig(f"{out_prefix}.read_len_distribution.boxplot.pdf", format='pdf', dpi=300)
    plt.savefig(f"{out_prefix}.read_len_distribution.boxplot.svg", format='svg')
    plt.close()




def align2bed(ReadsMapDF, outfile):
        df = ReadsMapDF.copy()
        df = df[df['read_type'].isin(['intra','inter'])]
        with open(outfile, 'w') as f:
                head_info = "\t".join(["#reference", "ref_start", "ref_end", "order"])
                f.write(f"{head_info}\n")
                for _, row in df.iterrows():
                        align_dict = row['align']
                        order = row['order'] 
                        for ref_id , align_list in align_dict.items():
                                for align in align_list:
                                        bed_info = "\t".join([ref_id, str(align[0]), str(align[1]), str(order)])
                                        f.write(f"{bed_info}\n")





def main():
        args = parse_args()

        ReadsMapDF = MappingStat(args.frag_bam, args.Qthreshold)
        total_raw_readnum, total_raw_baseslen =raw_data_info(args.seqkit_stat)
        output_read_type_stats(ReadsMapDF, total_raw_readnum, total_raw_baseslen, args.ID)
    
        order_read_type_read(ReadsMapDF,total_raw_readnum, f'{args.ID}.order_read_type_read.stat')
        order_read_type_base(ReadsMapDF,total_raw_baseslen, f'{args.ID}.order_read_type_base.stat')
        order_ref_num_read(ReadsMapDF, f'{args.ID}.order_ref_num_read.stat')
        order_ref_num_base(ReadsMapDF, f'{args.ID}.order_ref_num_base.stat')

        cutter_gap(ReadsMapDF, f'{args.ID}.gap.stat')
        align2bed(ReadsMapDF, f'{args.ID}.align.bed')
        #plot_read_len_distribution(ReadsMapDF, args.ID)
        #map_length(ReadsMapDF, f'{args.ID}.map_length.stat')



if __name__ == "__main__":
        main()
