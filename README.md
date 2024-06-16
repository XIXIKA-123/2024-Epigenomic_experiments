# 2024学年表观组学实验
## 简介
该项目基于2024学年表观组学实验，通过分析表观组学数据（基于peak.bed数据），获得Promoter区域以及Enhancer区域的关联基因富集结果，以及motif scan分析。该实验涉及R语言、GREAT(http://great.stanford.edu/public/html/index.php)的使用、IGV可视化，以及Motif scan的两种方法（Monalisa、Motifscan）。
## 实验一
- 通过编程获取两个peak.bed文件的交集部分，并且返回来自于peakA.bed文件中的peak信息，可参考samtools中的相关功能。
- 代码见lab1.R
## 实验二
- 鉴定Promoter区域和Enhancer区域，鉴定规则基于下述假设：获取人类的TSS位点信息，以TSS位点为中心，向上下游各延伸2kb设置一个windows，若该windows与实验1中获取的交集peak存在交集，则认为该区域为Promoter区域，否则为Enhancer区域。
- 代码见lab2.R
## 实验三
- 对Promoter区域及Enhancer区域peak对应的基因进行富集分析。
- 该实验使用GREAT(http://great.stanford.edu/public/html/index.php)进行，参考基因组版本为hg38。富集可视化代码见lab3.R
## 实验四
- Motifscan可以获取转录因子以及对应motif序列和信息，该部分使用Monalisa和Motifscan两种方法进行。
