#!/usr/bin/env Rscript
# 加载所需的库
library(ggplot2)

# 从命令行参数获取输入文件名和输出文件名
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("请提供输入文件名和输出文件名。\nusage: Rscript order_type.r <input_file> <output_file>", call. = FALSE)
}
input_file <- args[1]
output_file <- args[2]

# 读取数据
file <- read.table(input_file,, header = TRUE, sep = "\t")
type <- factor(c(rep("intra", 9), rep("inter", 9)), levels = c("intra", "inter"))
read_order <- factor(rep(c(2, 3, 4, 5, 6, 7, 8, 9, ">=10"), 2), levels = c(">=10", 9, 8, 7, 6, 5, 4, 3, 2))
percent <- c(file$intra_read_num, file$inter_read_num)
data <- data.frame(type, read_order, percent)

# 创建图形
ggplot(data, aes(fill=read_order, y=type, x=percent)) + 
  geom_bar(position="stack", stat="identity",color="grey10",linewidth =0.6)+
  scale_fill_brewer(palette = "YlGn")+
  theme(text=element_text(size= 45 ),plot.title = element_text(size = 45))

# 保存图形为文件，关闭尺寸限制
ggsave(output_file, plot = last_plot(), device = "png", width = 25, height = 8, units = "in", limitsize = FALSE)
