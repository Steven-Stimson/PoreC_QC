#!/usr/bin/env python3
import bisect
import sys

def count_cut_sites(bed_file, cut_pos_file, n=4):
    # 读取cut.pos文件并按refid分组，排序pos
    cuts = {}
    with open(cut_pos_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts[0].startswith("#"):
                continue  # 跳过格式错误的行
            refid, pos = parts
            pos = int(pos)
            if refid not in cuts:
                cuts[refid] = []
            cuts[refid].append(pos)
    
    # 对每个refid的pos列表排序
    for refid in cuts:
        cuts[refid].sort()
    
    # 处理BED文件，计算每个区间的统计信息
    results = []
    with open(bed_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if parts[0].startswith("#"):
                continue  # 跳过格式错误的行
            refid, start_str, end_str = parts[0:3]
            start = int(start_str)
            end = int(end_str)
            
            # 初始化计数器
            count_internal = 0  # 窗口内的位点数（排除误差范围）
            endpoints_count = 0  # 起终点覆盖计数
            
            if refid in cuts:
                pos_list = cuts[refid]
                
                # ------ 修正窗口内的计数逻辑 ------
                # 新起始：start + n（排除起点误差范围）
                # 新结束：end - n（排除终点误差范围）
                # 确保窗口内部区域不与误差范围重叠
                new_start = start + n
                new_end = end - n
                
                # 如果新区域无效（start >= end），计数为0
                if new_start < new_end:
                    left_internal = bisect.bisect_right(pos_list, new_start)
                    right_internal = bisect.bisect_left(pos_list, new_end)
                    count_internal = max(0, right_internal - left_internal)
                else:
                    count_internal = 0
                
                # 检查起点是否被覆盖（start±n）
                start_low = start - n
                start_high = start + n
                left_start = bisect.bisect_left(pos_list, start_low)
                right_start = bisect.bisect_right(pos_list, start_high)
                if right_start > left_start:
                    endpoints_count += 1
                
                # 检查终点是否被覆盖（end±n）
                end_low = end - n
                end_high = end + n
                left_end = bisect.bisect_left(pos_list, end_low)
                right_end = bisect.bisect_right(pos_list, end_high)
                if right_end > left_end:
                    endpoints_count += 1
            
            # 组合结果（格式：refid start end endpoints_count internal_count）
            result_line = f"{line.strip()}\t{endpoints_count}\t{count_internal}"
            results.append(result_line)
    
    return results

# 使用示例：
bed_file = sys.argv[1]    # 替换为你的BED文件路径
cut_pos_file = sys.argv[2]  # 替换为你的cut.pos文件路径

output = count_cut_sites(bed_file, cut_pos_file, n=int(sys.argv[3]))  # 默认n=4

# 打印结果或保存到文件
for line in output:
    print(line)
