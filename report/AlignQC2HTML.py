#!/usr/bin/env python
import argparse
from bs4 import BeautifulSoup
import datetime

def replace_tables(soup, table_paths):
    """ 按顺序替换所有表格内容 """
    table_containers = soup.find_all('div', {'class': 'table-container'})
    for table_path, container in zip(table_paths, table_containers):
        with open(table_path, 'r', encoding='utf-8') as f:
            lines = f.read().splitlines()

        # 创建表格结构
        table = soup.new_tag('table')
        thead = soup.new_tag('thead')
        tbody = soup.new_tag('tbody')

        # 解析表头
        headers = lines[0].split('\t')
        tr = soup.new_tag('tr')
        for h in headers:
            th = soup.new_tag('th', string=h)
            tr.append(th)
        thead.append(tr)

        # 解析数据行
        for line in lines[1:]:
            parts = line.split('\t')
            tr = soup.new_tag('tr')
            for part in parts:
                td = soup.new_tag('td', string=part)
                tr.append(td)
            tbody.append(tr)

        table.append(thead)
        table.append(tbody)

        # 替换占位符
        container.find('table').replace_with(table)

def replace_images(soup, image_paths):
    """ 按顺序替换所有图片路径 """
    image_tags = soup.find_all('img')  # 找到所有<img>标签
    for img_tag, image_path in zip(image_tags, image_paths):
        img_tag['src'] = image_path  # 替换图片路径

def replace_entity_and_reference(soup, entity_id, reference, resolution):
    """ 替换样本名称和参考基因组 """
    # 替换样本名称
    sample_name_tag = soup.find('span', id="sample-name")
    if sample_name_tag:
        sample_name_tag.string = entity_id

    # 替换参考基因组
    reference_genome_tag = soup.find('span', id="reference-genome")
    if reference_genome_tag:
        reference_genome_tag.string = reference

    # 替换参考分辨率
    resolution_tag = soup.find('span', id="recommended-resolution")
    if resolution_tag:
        resolution_tag.string = resolution


def update_report_time(soup):
    """ 更新报告生成时间 """
    report_time_tag = soup.find('p', class_="lead")
    if report_time_tag:
        current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        report_time_tag.string = f"生成时间：{current_time}"



def main():
    parser = argparse.ArgumentParser(description='生成数据分析报告')
    parser.add_argument('-i', '--input', required=True, help='输入HTML模板文件路径')
    parser.add_argument('-o', '--output', required=True, help='输出HTML文件路径')
    parser.add_argument('--tables', nargs='+', required=True, help='表格文件路径列表')
    parser.add_argument('--images', nargs='+', required=True, help='图片文件路径列表')
    parser.add_argument('--EntityID', '-id', required=True, help='样本名称')
    parser.add_argument('--reference', '-ref', required=True, help='参考基因组')
    parser.add_argument('--resolution', '-res', required=True, help='分辨率')

    args = parser.parse_args()

    # 读取HTML模板
    with open(args.input, 'r', encoding='utf-8') as f:
        html_content = f.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # 替换表格内容
    replace_tables(soup, args.tables)

    # 替换图片路径
    replace_images(soup, args.images)

    # 替换样本名称和参考基因组
    replace_entity_and_reference(soup, args.EntityID, args.reference, args.resolution)

    # 更新报告生成时间
    update_report_time(soup)

    # 保存结果
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(str(soup))

if __name__ == '__main__':
    main()
