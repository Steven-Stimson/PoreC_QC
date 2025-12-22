#!/usr/bin/env python
import sys
import os
import base64
from bs4 import BeautifulSoup

def image2base64(html, image_dir):
    # 读取原始 HTML 文件
    with open(html, "r", encoding="utf-8") as f:
        html_content = f.read()

    soup = BeautifulSoup(html_content, 'html.parser')

    # 遍历所有图片文件
    for filename in os.listdir(image_dir):
        if filename.endswith((".svg", ".png", ".jpg")):
            file_path = os.path.join(image_dir, filename)
            # 转换为 Base64
            with open(file_path, "rb") as img_file:
                encoded_str = base64.b64encode(img_file.read()).decode('utf-8')
            # 根据扩展名生成正确的 data URL
            ext = filename.split('.')[-1].lower()
            if ext == 'svg':
                data_uri = f"data:image/svg+xml;base64,{encoded_str}"
            else:
                data_uri = f"data:image/{ext};base64,{encoded_str}"

            # 替换 HTML 中的路径
            old_src = f'src="{os.path.join(image_dir, filename)}"'
            new_src = f'src="{data_uri}"'
            html_content = html_content.replace(old_src, new_src)

    # 保存新 HTML（无需图片文件夹）
    with open(f'{html[:-5]}.embeded.html', "w", encoding="utf-8") as f:
        f.write(html_content)

def main():
    image2base64(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
