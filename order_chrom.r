#!/usr/bin/env Rscript
# 加载所需的库
library(ggplot2)

# 从命令行参数获取输入文件名和输出文件名
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  stop("请提供输入文件名和输出文件名。\nusage: Rscript order_chrom.r <input_file> <output_file>", call. = FALSE)
}
input_file <- args[1]
output_file <- args[2]

# 读取数据
file <- read.table(input_file,, header = TRUE, sep = "\t")
read_order<- factor(file$order,levels=c(2,3,4,5,6,7,8,9,">=10"))
chrom_num<- factor(file$chrom_num,levels=c(">=5",4,3,2,"intra"))
percent <- c(file$read_ratio)
sumdata <- data.frame(read_order,chrom_num,percent)

ggplot(sumdata, aes(fill=chrom_num, y=percent, x=read_order)) + 
  geom_bar(position="fill", stat="identity",color="grey10",linewidth=0.6)+
  scale_fill_brewer(palette = "YlGn")+
  theme(text=element_text(size= 45 ),plot.title = element_text(size = 45))

ggsave(output_file, plot = last_plot(), device = "png", width = 25, height = 21, units = "in", limitsize = FALSE)
