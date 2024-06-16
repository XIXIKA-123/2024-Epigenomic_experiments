# @Author:  Qian jing 
# @Email: qianjing03223@163.com

# 加载所需的R包
library("devtools") # 用于安装和加载其他包
library("monaLisa") # 用于motif富集分析和可视化
library("JASPAR2022") # 提供JASPAR数据库中的motif信息
library("rtracklayer") # 用于处理基因组注释和区域
library("TFBSTools") # 用于转录因子结合位点分析
library("SummarizedExperiment") # 用于存储和操作高通量的基因组数据
library("BSgenome.Hsapiens.UCSC.hg38") # 加载人类基因组参考序列


###---准备数据---###
# 读取bed文件中的peak数据
bed_file = "promoter.bed"
bed_peaks = import.bed(bed_file)

# 调整peak的大小，使它们中心对齐，并具有相同的宽度
single_set = resize(bed_peaks, 
                    width = median(width(bed_peaks)), # 使用peak宽度的中位数
                    fix = "center") # 以中心固定

# 从参考基因组中提取peak区域的序列
single_seqs = getSeq(BSgenome.Hsapiens.UCSC.hg38, single_set)

# 获取JASPAR2022数据库中的PWM矩阵，指定为脊椎动物
pwms = getMatrixSet(x = JASPAR2022,
                    opts = list(
                      matrixtype = "PWM", 
                      tax_group = "vertebrates" # 指定为脊椎动物的矩阵
                    ))

###---Motif富集分析---###
# 计算序列中的motif富集情况，并与基因组背景进行比较
sinse = calcBinnedMotifEnrR(seqs = single_seqs,
                            pwmL = pwms,
                            background = "genome", # 使用基因组作为背景
                            genome = BSgenome.Hsapiens.UCSC.hg38, # 指定参考基因组
                            genome.regions = NULL, 
                            genome.oversample = 2, # 从整个基因组中进行采样
                            BPPARAM = BiocParallel::SerialParam(RNGseed = 42), # 设置并行计算参数
                            verbose = TRUE) 

# 在peak区域寻找与PWM矩阵匹配的motif
res = findMotifHits(query = pwms,
                    subject = single_seqs,
                    min.score = 6.0, # 设置最小分数阈值
                    method = "matchPWM", # 使用匹配PWM的方法
                    BPPARAM = BiocParallel::SerialParam())

# 将富集分析的结果整理成数据框，并过滤出p值大于6的条目
motif_pvalue = data.frame(tf = rownames(assay(sinse, "negLog10Padj")),
                          p = assay(sinse, "negLog10Padj")[,1]) %>%
  na.omit() %>% # 移除含有NA值的行
  dplyr::filter(p > 6) %>% # 过滤p值大于6的条目
  dplyr::arrange(desc(p)) # 按p值降序排列

# 将motif的p值结果写入csv文件
write.csv(motif_pvalue, "promoter_motifs.csv")

###---可视化motif富集分析结果---###
# 选择前50个显著的motif进行热图可视化
plotMotifHeatmaps(x = sinse[motif_pvalue$tf[1:50],],
                  which.plots = c("log2enr", "negLog10Padj"), # 选择绘制的图形类型
                  width = 1.8, # 热图的宽度
                  maxEnr = 2, # 最大富集倍数
                  maxSig = 10, # 最大显著性水平
                  cluster = T, # 进行聚类
                  show_dendrogram = T, # 显示树状图
                  show_seqlogo = T) # 显示序列标志
