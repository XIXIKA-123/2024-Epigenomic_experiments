# @Author:  Qian jing 
# @Email: qianjing03223@163.com

setwd("D:/本科资料/大三下/表观组学/实验1")

peakOverlap <- function(file_A, file_B){
  # bed文件没有表头，因此header为F
  A <- read.table(file_A, header = FALSE, col.names = c("chr", "start", "end"))
  B <- read.table(file_B, header = FALSE, col.names = c("chr", "start", "end"))
  
  # 检索策略:按照染色体条目进行分类,返回A的peak，共四种情况
  overlap <- A[sapply(1:nrow(A), function(i){
    chr_A <- A$chr[i]
    start_A <- A$start[i]
    end_A <- A$end[i]
    idx <- which(B$chr == chr_A & (
      (start_A <= B$start & end_A >= B$end) | # B完全在A内
        (start_A >= B$start & start_A <= B$end) | # A的start在B内
        (end_A >= B$start & end_A <= B$end) | # A的end在B内
        (B$start <= start_A & B$end >= end_A) # A完全在B内
    ))
    length(idx) > 0
  }), ]
  
  # 保存为bed文件，分隔符为制表符
  write.table(overlap, "res.bed", sep = '\t',row.names = FALSE, col.names = FALSE, quote = FALSE)
}
# 运行函数，获得res.bed，文件中为来自于peakA.bed文件的交集peak信息
peakOverlap(file_A = "737_10_5.bed",file_B = "sg_10_5.bed")

