# @Author:  Qian jing 
# @Email: qianjing03223@163.com

setwd("D:/本科资料/大三下/QJ/表观组学/实验3")

GO_BP <- read.table("GOBP.tsv",sep="\t",header=T)
GO_CC <- read.table("GOCC.tsv",sep="\t",header=T)
GO_MF <- read.table("GOMF.tsv",sep="\t",header=T)

# 取前十个term
GO_BP <- GO_BP[1:10,]
GO_CC <- GO_CC[1:10,]
GO_MF <- GO_MF[1:10,]

# 添加-log10pvalue列
GO_BP$log10pvalue <- -log10(GO_BP$BinomBonfP)
GO_CC$log10pvalue <- -log10(GO_CC$BinomBonfP)
GO_MF$log10pvalue <- -log10(GO_MF$BinomBonfP)

# 合并数据并添加分类标签
GO_BP$Category <- "BP"
GO_CC$Category <- "CC"
GO_MF$Category <- "MF"
GO_combined <- rbind(GO_BP, GO_CC, GO_MF)
#####
# 绘制柱状图
# 设置颜色
colors <- c("BP" = "#ea8379", "CC" = "#7daee0", "MF" = "#b395bd")
bar_colors <- colors[GO_combined$Category]
# 打开PDF设备
pdf("GO_Term_Enrichment.pdf", width = 12, height = 7)

# 设置图形参数，绘制柱状图
par(mar = c(5, 12, 4, 10))  # 增加左边和右边的边距以放置图例和标签
p <- barplot(GO_combined$log10pvalue,
             names.arg = GO_combined$Desc,
             horiz = TRUE,
             col = bar_colors,
             border = "white",
             xlab = "-log10(padj)",
             main = "GO Term Enrichment",
             las = 1,  # 将纵坐标标签水平显示
             cex.names = 0.6,  # 调整标签字体大小，使其更小
             xlim = c(0, 50))  # 设置x轴范围为0到50

# 添加图例到图表外部，设置字体加粗和调整大小
legend("topright", 
       legend = names(colors), fill = colors, border = "white", 
       cex = 0.8, text.font = 2, xpd = TRUE, box.lwd = 2)

# 关闭PDF设备
dev.off()
