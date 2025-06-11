#!/usr/bin/env python3
import pysam
from collections import Counter
import argparse 
from collections import defaultdict
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
                    prog='inter&intra reads statistics',
                    description='get intra&inter&align_one reads informations')

parser.add_argument("--out_prefix", type=str, required=True, help="Prefix for the output files")
parser.add_argument("--bam", type=str, nargs="+", required=True, help="List of BAM files")
parser.add_argument("--binsize", type=int, default=50000, help="Binsize or resolution")
parser.add_argument("--mapq_threshold", type=int, default=10, help="MAPQ threshold for filtering")

args = parser.parse_args()

out_prefix = args.out_prefix
filepath = args.bam
binsize = args.binsize
mapq_threshold = args.mapq_threshold

#binsize = 50000
#mapq_threshold = 10

def read_type_classify(value,mapq_th=0):
    if mapq_th==0:
        read_len = np.sum(value["frag_len"])
        unique_chrom = set(x for x in value['chrom'] if x != '*')
        if len(unique_chrom)==0: read_tpye = "unmap" ##这是没有比对上的reads
        if len(unique_chrom)==1:
            if sum(1 for x in value['chrom'] if x != '*') == 1:
                read_tpye = "align_one"
            elif sum(1 for x in value['chrom'] if x != '*') > 1:
                read_tpye = "intra"
        if len(unique_chrom)>1: read_tpye = "inter" 
        return(read_len,read_tpye)
    else:
        valid_indices = []
        for idx, mapq in enumerate(value['mapq']):
            if mapq != '*' and int(mapq) >= mapq_th:
                valid_indices.append(idx)
        filtered_frag_len = [value['frag_len'][idx] for idx in valid_indices]
        filtered_chrom = [value['chrom'][idx] for idx in valid_indices]
        #filtered_bin = [value['bin'][idx] for idx in valid_indices]
        read_len = np.sum(filtered_frag_len)
        unique_chrom = set(filtered_chrom)
        if len(filtered_chrom)==0:read_tpye = "unmap"
        elif len(filtered_chrom)==1:
            read_tpye = "align_one"
        elif len(filtered_chrom) != 1 and len(unique_chrom)==1:
                read_tpye = "intra"
        elif len(filtered_chrom) != 1 and len(unique_chrom)>1:
            read_tpye = "inter" 
        return(read_len,read_tpye)

def count_prefixes(value):
    prefix_list = []
    for a in value:
        prefix = a.split('_')[0]
        prefix_list.append(prefix)
    counts = (Counter(prefix_list)).most_common()  ##返回的是一个[('2', 3), ('11', 1), ('4', 1), ('3', 1)]list里面是元组的结构,元组第一个是rst的编号，第二个是该rst的bin数
    return counts

def bin_set_stat(value, mapq_th=0):
    if mapq_th == 0:
        unique_bin = set(x for x in value['bin'] if x != '*')
        if unique_bin:  # 检查 unique_bin 是否非空
            read_whloe_bincount = len(unique_bin)
            inonechrom_maxcount = count_prefixes(unique_bin)[0][1]
            return read_whloe_bincount, inonechrom_maxcount
        else:
            return 0, 0
    else:
        valid_indices = [idx for idx, mapq in enumerate(value['mapq']) if mapq != '*' and int(mapq) >= mapq_th]
        if valid_indices:
            filtered_bin = [value['bin'][idx] for idx in valid_indices]
            unique_bin = set(filtered_bin)
            if unique_bin:
                read_whloe_bincount = len(unique_bin)
                inonechrom_maxcount = count_prefixes(unique_bin)[0][1]
                return read_whloe_bincount, inonechrom_maxcount
            else:
                return 0, 0
        else:
            return 0, 0 

def binset_up5merge(count_list):
    more5_count_list = [item for item in count_list if item[0] >= 5]
    keep_count_list = [item for item in count_list if item[0] < 5]
    more5_count = sum([item[1] for item in more5_count_list])
    keep_count_list.append((5,more5_count))
    return keep_count_list

def order_chrom_get(value,mapq_th):
    valid_indices = []
    order = sum(1 for x in value['mapq'] if x != '*' and x>=mapq_threshold)
    for idx, mapq in enumerate(value['mapq']):
        if mapq != '*' and int(mapq) >= mapq_th:
            valid_indices.append(idx)
    filtered_chrom_num = len(set([value['chrom'][idx] for idx in valid_indices]))
    return order , filtered_chrom_num

read_type_results = []
original_binset_results = []
filter_binset_results = []
order_type_results = []
order_chrom_results = []

for file in filepath:
    info_dict = defaultdict(lambda: {'mapq':[], 'frag_len': [], 'chrom': [] , 'bin':[]})
    with pysam.AlignmentFile(file, 'rb') as openfile:
        references = openfile.references
        for read in openfile:
            read_name = (read.query_name).split(':')[0]
            if read.flag== 4:
                mapq = "*"
                frag_len = read.query_length
                chrom_name = "*"
                bin_ID = "*"
                info_dict[read_name]['mapq'].append(mapq)
                info_dict[read_name]['frag_len'].append(frag_len)
                info_dict[read_name]['chrom'].append(chrom_name)
                info_dict[read_name]['bin'].append(bin_ID)
            else:
                mapq = read.mapq
                frag_len = read.query_length
                chrom_name = references[read.reference_id][3:] ##取染色体的编号，去掉前面的Chr
                median = (read.reference_start+read.reference_end)/2
                bin_ID = '{}_{}'.format(chrom_name, int(np.ceil(median/ binsize))) #ceil向上取整
                info_dict[read_name]['mapq'].append(mapq)
                info_dict[read_name]['frag_len'].append(frag_len)
                info_dict[read_name]['chrom'].append(chrom_name)
                info_dict[read_name]['bin'].append(bin_ID)
    read_type_dict = {"unmap":[],"align_one":[],"intra":[],"inter":[]}
    filter_read_type_dict = {"unmap":[],"align_one":[],"intra":[],"inter":[]}
    inwholechrom_bincount = [] ##仅仅获得的是reads的frag对应的去冗余的bin数
    inonechrom_maxbincount = []
    filter_inwholechrom_bincount = []
    filter_inonechrom_maxbincount = []
    limit=10
    limit_str='>=%d' % limit
    rangelist=list(range(2,limit))
    rangelist.append(limit_str)
    #intra_order_type={limit_str:[]}
    #inter_order_type={limit_str:[]}
    order_type = {
        "intra":{
            2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],limit_str:[]
        },
        "inter":{
            2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],limit_str:[]
        }
    }
    for read_name, value in info_dict.items():
        #read_type_classify过滤前
        read_len,read_tpye = read_type_classify(value)
        read_type_dict[read_tpye].append(read_len)
        #read_type_classify过滤后
        filter_read_len,filter_read_tpye = read_type_classify(value,mapq_threshold)
        filter_read_type_dict[filter_read_tpye].append(filter_read_len)        
        #bincount过滤前
        read_whloe_bincount,inonechrom_maxcount = bin_set_stat(value)
        inwholechrom_bincount.append(read_whloe_bincount)
        inonechrom_maxbincount.append(inonechrom_maxcount)
        #bincount过滤后
        filter_read_whloe_bincount,filter_inonechrom_maxcount = bin_set_stat(value,mapq_threshold)
        filter_inwholechrom_bincount.append(filter_read_whloe_bincount)
        filter_inonechrom_maxbincount.append(filter_inonechrom_maxcount)
        ##过滤后堆叠图
        if filter_read_tpye=="intra" or filter_read_tpye=="inter":
            order , filtered_chrom_num = order_chrom_get(value,mapq_threshold)
            if order < limit : order_type[filter_read_tpye][order].append(filtered_chrom_num)
            elif order >= limit : order_type[filter_read_tpye][limit_str].append(filtered_chrom_num)   
    ##reads_type输出汇总
    align_one_readnum=len(read_type_dict["align_one"])
    intra_readnum=len(read_type_dict["intra"])
    inter_readnum=len(read_type_dict["inter"])
    mapping_readnum = align_one_readnum + intra_readnum + inter_readnum
    align_one_basenum=sum(read_type_dict["align_one"])
    intra_basenum=sum(read_type_dict["intra"])
    inter_basenum=sum(read_type_dict["inter"])   
    mapping_basenum = align_one_basenum + intra_basenum + inter_basenum
    filter_align_one_readnum=len(filter_read_type_dict["align_one"])
    filter_intra_readnum=len(filter_read_type_dict["intra"])
    filter_inter_readnum=len(filter_read_type_dict["inter"])
    filter_mapping_readnum = filter_align_one_readnum + filter_intra_readnum + filter_inter_readnum
    filter_align_one_basenum=sum(filter_read_type_dict["align_one"])
    filter_intra_basenum=sum(filter_read_type_dict["intra"])
    filter_inter_basenum=sum(filter_read_type_dict["inter"])   
    filter_mapping_basenum = filter_align_one_basenum + filter_intra_basenum + filter_inter_basenum
    read_type_result = {
        "file_name": file,
        "align_one_read": align_one_readnum,
        "intra_read": intra_readnum,
        "inter_read": inter_readnum,
        "mapping_read": mapping_readnum,
        "align_one_base": align_one_basenum,
        "intra_base": intra_basenum,
        "inter_base": inter_basenum,
        "mapping_base": mapping_basenum,
        "filter_align_one_read": filter_align_one_readnum,
        "filter_intra_read": filter_intra_readnum,
        "filter_inter_read": filter_inter_readnum,
        "filter_mapping_read": filter_mapping_readnum,
        "filter_align_one_base": filter_align_one_basenum,
        "filter_intra_base": filter_intra_basenum,
        "filter_inter_base": filter_inter_basenum,
        "filter_mapping_base": filter_mapping_basenum
    }
    read_type_results.append(read_type_result)
    ##reads_bin输出汇总
    a = sorted(binset_up5merge(Counter(inwholechrom_bincount).most_common())) #reads的frag对应的去冗余的bin数及其在reads中的频数 即bin数→频数
    b = sorted(binset_up5merge(Counter(inonechrom_maxbincount).most_common()))
    filter_a = sorted(binset_up5merge(Counter(filter_inwholechrom_bincount).most_common())) #reads的frag对应的去冗余的bin数及其在reads中的频数 即bin数→频数
    filter_b = sorted(binset_up5merge(Counter(filter_inonechrom_maxbincount).most_common()))
    original_binset_result={
        "file_name": file,
        "inwholechrom_bincount0":a[0][1],
        "inwholechrom_bincount1":a[1][1],
        "inwholechrom_bincount2":a[2][1],
        "inwholechrom_bincount3":a[3][1],
        "inwholechrom_bincount4":a[4][1],
        "inwholechrom_bincount>=5":a[5][1],
        "inonechrom_maxbincount0":b[0][1],
        "inonechrom_maxbincount1":b[1][1],     
        "inonechrom_maxbincount2":b[2][1],
        "inonechrom_maxbincount3":b[3][1],
        "inonechrom_maxbincount4":b[4][1],
        "inonechrom_maxbincount>=5":b[5][1]
    }
    original_binset_results.append(original_binset_result)
    filter_binset_result={
        "file_name": file,
        "inwholechrom_bincount0":filter_a[0][1],
        "inwholechrom_bincount1":filter_a[1][1],
        "inwholechrom_bincount2":filter_a[2][1],
        "inwholechrom_bincount3":filter_a[3][1],
        "inwholechrom_bincount4":filter_a[4][1],
        "inwholechrom_bincount>=5":filter_a[5][1],
        "inonechrom_maxbincount0":filter_b[0][1],
        "inonechrom_maxbincount1":filter_b[1][1],     
        "inonechrom_maxbincount2":filter_b[2][1],
        "inonechrom_maxbincount3":filter_b[3][1],
        "inonechrom_maxbincount4":filter_b[4][1],
        "inonechrom_maxbincount>=5":filter_b[5][1]
    }
    filter_binset_results.append(filter_binset_result)
    ##order_type 过滤后
    for order_ in rangelist:
        order_type_result={
            "file_name": file,
            "order":order_,
            "intra_read_num":len(order_type["intra"][order_]),
            "inter_read_num":len(order_type["inter"][order_])
        }
        order_type_results.append(order_type_result)
    ##order_chrom 过滤后
    pass_list=[2,3,4]
    for order_ in rangelist:
        order_total_num=len(order_type["intra"][order_])+len(order_type["inter"][order_])
        intra_a_info=dict(Counter(order_type["intra"][order_]))
        inter_a_info=dict(Counter(order_type["inter"][order_]))
        pass_inter_a_info={key: value for key, value in inter_a_info.items() if key not in pass_list}
        order_chrom_result = {
            "file_name": file,
            "order":order_,
            "chrom_num":"intra",
            "read_num":intra_a_info[1]
        }
        order_chrom_results.append(order_chrom_result)
        for aa in range(2,5):
            if aa in inter_a_info.keys():
                order_chrom_result = {
                    "file_name": file,
                    "order":order_,
                    "chrom_num":aa,
                    "read_num":inter_a_info[aa]
                }
                order_chrom_results.append(order_chrom_result)
            else:
                order_chrom_result = {
                    "file_name": file,
                    "order":order_,
                    "chrom_num":aa,
                    "read_num":0.00
                }
                order_chrom_results.append(order_chrom_result)
        up5_num=sum(pass_inter_a_info.values())
        order_chrom_result = {
            "file_name": file,
            "order":order_,
            "chrom_num":">=5",
            "read_num":up5_num
        }
        order_chrom_results.append(order_chrom_result)

outfile1 = '01.{}.read_type.txt'.format(out_prefix)
outfile2 = '02.{}.rawdata.bincount.txt'.format(out_prefix)
outfile3 = '03.{}.filterdata.bincount.txt'.format(out_prefix)
outfile4 = '04.{}.order_type.split.txt'.format(out_prefix)
outfile5 = '05.{}.order_chrom.split.txt'.format(out_prefix)
#outfile6 = '06.{}.order_type.merge.txt'.format(out_prefix)
#outfile9 = '07.{}.order_chrom.merge.txt'.format(out_prefix)

df1 = pd.DataFrame(read_type_results)
df2 = pd.DataFrame(original_binset_results)
df3 = pd.DataFrame(filter_binset_results)
df4 = pd.DataFrame(order_type_results)
df5 = pd.DataFrame(order_chrom_results)
df1.to_csv(outfile1, index=False,sep="\t")
df2.to_csv(outfile2, index=False,sep="\t")
df3.to_csv(outfile3, index=False,sep="\t")
df4.to_csv(outfile4, index=False,sep="\t")
df5.to_csv(outfile5, index=False,sep="\t")

'''
#merge order_type_results
df6 = df4.groupby("order").agg({'intra_read_num': 'sum', 'inter_read_num': 'sum'}).apply(lambda x: round(x *100 / sum(df1["filter_mapping_read"]),2))
#merge order_chrom_results
df7 = df5.groupby(["order","chrom_num"]).agg({'read_num': 'sum'}).reset_index()
df8 = df7.groupby("order").agg({'read_num': 'sum'})
df9 = pd.merge( df7,df8, on="order", how="outer")
df9['read_ratio'] = round(df9['read_num_x'] / df9['read_num_y']*100,2)
df6.to_csv(outfile6, sep="\t")
df9.to_csv(outfile9, index=False,sep="\t")
'''
