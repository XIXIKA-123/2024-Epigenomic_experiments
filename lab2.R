# @Author:  Qian jing 
# @Email: qianjing03223@163.com

setwd("D:/本科资料/大三下/表观组学/实验2")

# 染色体规模
chr_size <- read.table("hg38.chrom.sizes",header = T)
# 获取TSS位点
HM <- read.table("hg38.gene.txt",header=T)

# 根据正负链，即正链的start和负链的end作为TSS位点
plus_strand <- HM[HM$strand == "+", c("chrom", "txStart")]
colnames(plus_strand) <- c("chrom", "TSS")

minus_strand <- HM[HM$strand == "-", c("chrom", "txEnd")]
colnames(minus_strand) <- c("chrom", "TSS")

TSS_info <- rbind(plus_strand, minus_strand)
### 该步骤中即限制了chrom的size（start-2000如果小于0则取1，end+2000大于size则取size）
### 同时在match中过滤了非标准命名的chrom
TSS_extend <- data.frame(chrom = TSS_info$chrom,
                         TSS_start = pmax(TSS_info$TSS - 2000, 1),
                         TSS_end = pmin(TSS_info$TSS + 2000, 
                                        chr_size$size[match(TSS_info$chrom, chr_size$chr)]))
# 删除含有NA值的行
TSS_extend <- na.omit(TSS_extend)
TSS_extend <- unique(TSS_extend)
write.table(TSS_extend,"TSS_extend.bed",sep="\t",row.names = FALSE, col.names = FALSE, quote = FALSE)

# 实验1中的结果
lab1_bed <- read.table("res.bed")

peakNonOverlap <- function(file_A, file_B) {
  A <- read.table(file_A, header = FALSE, col.names = c("chr", "start", "end"))
  B <- read.table(file_B, header = FALSE, col.names = c("chr", "start", "end"))
  
  # 获取重叠区域
  overlap <- A[sapply(1:nrow(A), function(i) {
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
  
  # 写入重叠区域，这里对应的就是promoter区域
  write.table(overlap, "promoter.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 获取非重叠区域，即取非交集部分，这里对应的是enhancer区域
  nonOverlap <- A[!(row.names(A) %in% row.names(overlap)), ]
  write.table(nonOverlap, "enhancer.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

tss <- read.table("TSS_extend.bed")
peakNonOverlap(file_A = "res.bed",file_B = "TSS_extend.bed")

enhancer <- read.table("enhancer_edit.bed",header = F)
# 去除非标准的染色体，保存为enhancer_edit.bed
enhancer_edit <- enhancer[enhancer$V1 %in% tss$V1,]
write.table(enhancer_edit, "enhancer_edit.bed", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


