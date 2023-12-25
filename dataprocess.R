remove(list = ls()) #一键清空
options(stringsAsFactors = FALSE) 
Sys.setenv(LANGUAGE = "en") #设置显示语言
#设置工作目录
setwd("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\")

# 下载R包 --------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
BiocManager::install("pheatmap")#通过BiocManager安装我们常用的R包

# 整理TCGA下载的RNA-seq数据 --------------------------------------------------------------------
setwd("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\")
library("rjson")
json <- jsonlite::fromJSON("metadata.cart.json")
View(json)
#id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  

#获取gdc_download文件夹下的所有TSV表达文件的 路径+文件名，这里路径要改为自己的路径
#count_file <- list.files('gdc_download',pattern = '*.tsv',recursive = TRUE)  
count_file <- list.files('gdc_download',pattern = '*gene_counts.tsv',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

matrix = data.frame(matrix(nrow=60660,ncol=0))    #60660是tsv文件中的基因数，恒定是60660个基因，不用改代码。
for (i in 1:length(count_file_name)){
  path = paste0('gdc_download//',count_file[i])   #这里要改路径
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[2,]
  data <-data[-c(1:6),]
  data <- data[3] #取出unstranded列（得到COUNT矩阵），若想提取fpkm-unstranded则改为data[7]，fpkm-up-unstranded改为data[8]
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  matrix <- cbind(matrix,data)
}

# 设置Gene Symbol为列名的矩阵（前面得到的是Ensembl ID）
path = paste0('gdc_download//',count_file[1])
data<- as.matrix(read.delim(path,fill = TRUE,header = FALSE,row.names = 1))
gene_name <-data[-c(1:6),1]
gene_matrix <- cbind(gene_name,matrix)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果
gene_matrix <- aggregate( . ~ gene_name,data=gene_matrix, max)    
#将gene_name列设为行名
rownames(gene_matrix) <- gene_matrix[,1]
gene_matrix <- gene_matrix[,-1]
write.csv(gene_matrix,'Gene_counts.csv',row.names = TRUE)

# 处理RNA-seq数据，去重，分组 --------------------------------------------------------------------
library(openxlsx)
library(tidyverse)
library(limma)
library(readr)

#1.1 读取表达量数据，如果工作区数据清空就读取，没就接着用
setwd("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\")
filepath_gene <- "Gene_counts.csv"
gene_matrix <- read_csv(filepath_gene)
dim(gene_matrix) 
gene_matrix[6:8,1:3] 
gene_matrix1 <- gene_matrix     #留个备份

#1.2 利用limma包对重复的基因取平均值合并
gene_matrix <- as.data.frame(avereps(gene_matrix[,],ID = rownames(gene_matrix)))

#1.3 行名重命名：取行名的1-16位
colnames(gene_matrix) <- substring(colnames(gene_matrix),1,16) %>% gsub("\\.","-",.)
#colnames(gene_matrix) <- gsub("\\.","-",colnames(gene_matrix))
gene_matrix[1:4,1:4]

#1.4 去除缺失值，去除50%样本中没有表达的基因
# 检查是否有缺失值
if (any(is.na(gene_matrix))) {
  # 如果存在缺失值，可以选择去除缺失值或进行其他处理
  gene_matrix <- na.omit(gene_matrix)
}
#计算每行从第1列到最后一列为表达值为0的样本占比
zero_percentage <- rowMeans(gene_matrix[, 1:ncol(gene_matrix)] == 0)
# 设置阈值为0.5，筛选满足条件基因
gene_matrix <- gene_matrix[zero_percentage < 0.5, ]
##去除低表达的基因，也可以卡其它阈值
gene_matrix <- gene_matrix[rowMeans(gene_matrix)>=1, ] 

#1.5 去除vial为B的数据
# B：很少数，表示福尔马林固定石蜡包埋组织，测序分析的效果不佳，删掉01B的样本数据
table(substring(colnames(gene_matrix),16))
Vial_A <- colnames(gene_matrix)[as.character(substr(colnames(gene_matrix),16,16)) == "A"]
gene_matrix <- gene_matrix[,Vial_A]
colnames(gene_matrix) <- substr(colnames(gene_matrix), 1, 15)
gene_matrix[6:8,1:3] 

#1.6 去除批次效应（有待完成）


'"#1.7 根据表达矩阵病人的ID进行分组
TCGA_group_list <- ifelse(as.numeric(substring(colnames(gene_matrix),14,15)) < 10,
                          "Tumor","Normal") %>% 
  factor(.,levels = c("Normal","Tumor"))
table(TCGA_group_list) # Normal 4   Tumor 8    "'  

#1.7 导出清洗好的数据
gene_matrix <- t(gene_matrix)
write.csv(gene_matrix,"C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_final.csv",quote = F,row.names = T)
write.table(gene_matrix, 'C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\Immune\\DATA.txt', sep="\t", quote =FALSE) #后面免疫浸润要用

# 识别mRNA，lncRNA和miRNA -----------------------------------------------------

# 2.1 读入EXCEL文件的gene注释信息，这里用到从biotype上下载好的EXCEL文件
mRNA_info <- read.xlsx("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA-seq\\Gene_info.xlsx",sheet = "mRNA_info")
lncRNA_info <- read.xlsx("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA-seq\\Gene_info.xlsx",sheet = "lncRNA_info")
miRNA_info <- read.xlsx("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA-seq\\Gene_info.xlsx",sheet = "miRNA_info")

#2.2 根据基因的注释信息，提取对应的表达矩阵
mRNA_matrix <- t(gene_matrix[rownames(gene_matrix) %in% mRNA_info$gene_name,])
dim(mRNA_matrix) 
write.csv(mRNA_matrix,"D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_mRNA.csv",quote = F,row.names = T)

lncRNA_matrix <- t(gene_matrix[rownames(gene_matrix) %in% lncRNA_info$gene_name,])
dim(lncRNA_matrix) 
write.csv(lncRNA_matrix,"D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_lncRNA.csv",quote = F,row.names = T)

miRNA_matrix <- t(gene_matrix[rownames(gene_matrix) %in% miRNA_info$gene_name,])
dim(miRNA_matrix) 
write.csv(miRNA_matrix,"D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_miRNA.csv",quote = F,row.names = T)

# 临床数据处理：划分ALN转移组和ALN非转移组 ----------------------------------------------------------------
setwd("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\")

#整理临床数据
library("rjson")
json <- jsonlite::fromJSON("RNA_seq\\metadata.cart.json")
View(json)
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))

clinical <- read.delim('Clinical_data\\clinical\\clinical.tsv',header = T)
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])
clinical_matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
clinical_matrix <- clinical_matrix[,-1]
write.csv(clinical_matrix,'Clinical_data\\clinical_matrix.csv',row.names = TRUE)


#分组：根据ALN转移与否分组：ALN_clinical和Not_ALN_clinical
library(dplyr)     #数据分析包，类似于python中的pandas

ALN_info <- c("N1a", "N1c", "N2a", "N3", "N3a", "N3b")
Not_ALN_info <- c("N0", "N0(i-)","N0 (i-)", "N0 (i+)","N0(i+)","N0 (mol+)", "N0(mol+)", "N1b", "N2b")

ALN_clinical <- clinical_matrix[clinical_matrix$ajcc_pathologic_n %in% ALN_info, ]
Not_ALN_clinical <- clinical_matrix[clinical_matrix$ajcc_pathologic_n %in% Not_ALN_info, ]

write.csv(ALN_clinical,"Clinical_data\\ALN_clinical.csv",quote = F,row.names = T)
write.csv(Not_ALN_clinical,"Clinical_data\\Not_ALN_clinical.csv",quote = F,row.names = T)

# 基因数据分组：②肿瘤组织；②③划分ALN组--------------------------------
library(readr)
library(dplyr)
filepath_gene <- "C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_mRNA.csv"
gene_matrix <- read_csv(filepath_gene)
gene_matrix <- as.data.frame(gene_matrix)
row.names(gene_matrix) <- gene_matrix[,1]
gene_matrix <- gene_matrix[, -1]
gene_matrix[6:8,1:3] 

##1.提取出肿瘤组织
#TCGA Barcode的第14和15位代表了样本的类型，01是实体瘤，11是正常样本
table(substring(colnames(gene_matrix),14,15))
library(stringr)
tumor <- colnames(gene_matrix)[as.integer(substr(colnames(gene_matrix),14,15)) == 01]
normal <- colnames(gene_matrix)[as.integer(substr(colnames(gene_matrix),14,15)) == 11]

#将tumor样本和正常样本按顺序储存到一个矩阵中
tumor_sample <- gene_matrix[,tumor]
normal_sample <- gene_matrix[,normal]

write.csv(tumor_sample,"C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_tumor.csv",quote = F,row.names = T)
#exprSet_by_group <- cbind(tumor_sample,normal_sample)
#group_list <- c(rep('tumor',ncol(tumor_sample)),rep('normal',ncol(normal_sample)))

##2.划分ALN组和Not_ALN组
file_path1 <- "C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\Clinical_data\\ALN_clinical.csv"
file_path2 <- "C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\Clinical-data\\Not_ALN_clinical.csv"

ALN_clinical <- read.csv(file_path1)
Not_ALN_clinical <- read.csv(file_path2)

library(dplyr)
colnames(tumor_sample) <- substring(colnames(tumor_sample),1,12)

ALN_ID <- intersect(colnames(tumor_sample), ALN_clinical$case_submitter_id)
Not_ALN_ID <- intersect(colnames(tumor_sample), Not_ALN_clinical$case_submitter_id)

ALN_tumor <- tumor_sample[, ALN_ID]
Not_ALN_tumor <- tumor_sample[, Not_ALN_ID]

write.csv(ALN_tumor,"RNA_seq\\ALN_tumor.csv",quote = F,row.names = T)
write.csv(Not_ALN_tumor,"RNA_seq\\Not_ALN_tumor.csv",quote = F,row.names = T)

#将ALN样本和Not_ALN样本按顺序储存到一个矩阵中
exprSet_by_group <- cbind(ALN_tumor,Not_ALN_tumor)
group_list <- c(rep('ALN',ncol(ALN_tumor)),rep('Not_ALN',ncol(Not_ALN_tumor)))

# 使用DESeq2进行差异分析 ----------------------------------------------------------
#DESeq2包
install.packages("DESeq2") #没有DESeq2包的可以先进行安装
library(DESeq2)
#分组
condition = factor(group_list)
coldata <- data.frame(row.names = colnames(exprSet_by_group), condition)
write.csv(coldata,"RNA_seq\\ALN_info.csv",quote = F,row.names = T)
dds <- DESeqDataSetFromMatrix(countData = round(exprSet_by_group),
                              colData = coldata,
                              design = ~condition)
head(dds)
dds$condition<- relevel(dds$condition, ref = "Not_ALN") # 指定哪一组作为对照组

#差异分析
dds <- DESeq(dds)
allDEG2 <- as.data.frame(results(dds))
save(allDEG2,file = 'RNA_seq\\countsallDEG2.Rdata')

#筛选显著性差异的基因
library(dplyr)
#这里使用log2FoldChange > 2 且adj.P.Val < 0.05的作为差异基因，可以根据需求改变阈值大小
#nrDEG_DESeq2_signif <- allDEG2 %>% filter(log2FoldChange > 1) %>% filter(padj < 0.05)
nrDEG_DESeq2_signif <- allDEG2 %>% filter(log2FoldChange > 1) %>% filter(padj < 0.1)
save(nrDEG_DESeq2_signif,file = 'RNA_seq\\nrDEG_DESeq2_signif.Rdata')


# 使用edgeR包进行差异分析  ---------------------------------------------------------

#install.packages("edgeR") #没有edgeR包的可以先进行安装
library(edgeR)
group_list = factor(group_list)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(exprSet_by_group)
colnames(design) <- levels(group_list)

##差异分析
DGElist <- DGEList(counts = exprSet_by_group, group = group_list)
## 使用cpm值对低表达量的基因进行过滤
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 前面做过过滤，这里可做，也可以不做
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)

fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)

#筛选显著性差异的基因
library(dplyr)
#这里使用logFC > 2 且FDR < 0.05的作为差异基因，可以根据需求改变阈值大小
nrDEG_edgeR_signif <- nrDEG_edgeR %>% filter(logFC > 2) %>% filter(FDR < 0.05)
save(nrDEG_edgeR_signif,file = 'nrDEG_edgeR_signif.Rdata')
write.csv(nrDEG_edgeR_signif,'D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\nrDEG_edgeR_signif.csv') 

# GO&KEGG分析 ---------------------------------------------------------------

library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(Rgraphviz)
library(pathview)

##读取差异表达基因，将基因ID从GENE_SYMBOL转换为ENTREZ_ID
#已有OrgDb的常见物种#
BiocManager::install("org.Hs.eg.db") 
library(org.Hs.eg.db)
#基因ID转换#
keytypes(org.Hs.eg.db) #查看所有可转化类型
genneid <- as.vector(rownames(nrDEG_edgeR_signif))
nrDEG_edgeR_signif$gene_name <- genneid
entrezid_all = mapIds(x = org.Hs.eg.db,  #id转换的比对基因组（背景基因）的物种，以人为例
                      keys = genneid, #将输入的gene_name列进行数据转换
                      keytype = "SYMBOL", #输入数据的类型
                      column = "ENTREZID")#输出数据的类型
entrezid_all  = na.omit(entrezid_all)  #na省略entrezid_all中不是一一对应的数据情况
entrezid_all = data.frame(entrezid_all) #将entrezid_all变成数据框格式
head(entrezid_all)

###GO富集分析###
GO_enrich <- enrichGO(gene = entrezid_all[,1], #表示前景基因，即待富集的基因列表;[,1]表示对entrezid_all数据集的第1列进行处理
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", #设定读取的gene ID类型
                      ont = "ALL", #可以输入CC/MF/BP/ALL，ALL包括 Biological Process,Cellular Component,Mollecular Function三部分
                      #universe = 背景数据集 # 表示背景基因，无参的物种选择组装出来的全部unigenes作为背景基因；有参背景基因则不需要。
                      pvalueCutoff = 0.1,#设定p值阈值
                      qvalueCutoff = 0.1, #设定q值阈值，阈值设置太严格可导致筛选不到基因。可指定 1 以输出全部
                      pAdjustMethod = "BH",
                      readable = T) #是否将基因ID映射到基因名称。
go_enrich <- data.frame(GO_enrich) #将GO_enrich导成数据框格式
## 自带的绘图
barplot(GO_enrich,showCategory = 20)
dotplot(GO_enrich, showCategory=15)

#数据导出#
write.csv(go_enrich,'C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\GO_KEGG\\GO_enrich.csv') 

###GO/KEGG富集结果可视化###
#数据载入与处理#
library(ggplot2)
GO_result <- go_enrich
# GO三种类别，每种选择显著性最高的12个展示出来
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -p.adjust)

# 纵向柱状图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="Top 10 Enriched GO Terms") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()

#横向柱状图#
ggplot(go_enrichment_pathway, 
       aes(x=reorder(Description, -Count),y=Count, fill=ONTOLOGY)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen")) + #柱状图填充颜色
  facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  labs(x="GO Term", y="Gene_Number", title="Top 10 Enriched GO Terms")+
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

# 绘制气泡图
ggplot(go_enrichment_pathway, aes(x=reorder(Description, Count), y=Count)) +
  geom_point(aes(size=Count,color=-log10(p.adjust))) +
  scale_size_continuous(range=c(1, 10)) +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  theme_minimal() +
  scale_color_gradient(low = "pink",high ="red")+
  labs(color=expression(-log10(p.adjust),size="Count"), 
       x="Gene Ratio",y="Gene_Number",title="GO Enrichment")+
  theme_bw()


###--------KEGG富集分析--------###
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], #即待富集的基因列表
                         keyType = "kegg",
                         pAdjustMethod = 'BH',  #指定p值校正方法
                         organism= "human",  #hsa，可根据你自己要研究的物种更改，可在https://www.kegg.jp/brite/br08611中寻找
                         qvalueCutoff = 0.32, #指定 p 值阈值（可指定 1 以输出全部）
                         pvalueCutoff=0.32) #指定 q 值阈值（可指定 1 以输出全部）
kegg_enrich <- data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\GO_KEGG\\KEGG_enrich.csv') #数据导出
## 自带的绘图
barplot(KEGG_enrich, drop = TRUE, showCategory = 15,color = "p.adjust",title = "KEGG Pathway")
## 气泡图
dotplot(KEGG_enrich, showCategory=15)

###GO/KEGG富集结果可视化###
# 根据p-value值排序，选择显著的通路或功能（这里设定阈值为0.05）
significant_pathways <- subset(kegg_enrich, p.adjust < 0.5)

# 绘制基于ggplot2的条形图
ggplot(significant_pathways, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "KEGG Pathway Enrichment Analysis",
       x = "Pathway",
       y = "-log10(P-value)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

# 设置颜色的阶梯和对应的颜色
significant_pathways$pvalue_group <- cut(significant_pathways$p.adjust, breaks = c(0, 0.001, 0.01, 0.05), labels = c("p < 0.001", "0.001 <= p < 0.01", "0.01 <= p < 0.05"))

# 绘制基于ggplot2的气泡图
ggplot(significant_pathways, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
  geom_point(aes(size = Count, color = pvalue_group), alpha = 0.7) +
  scale_size_continuous(range = c(1, 10)) +
  labs(title = "KEGG Pathway Enrichment Analysis",x = "Pathway",y = "-log10(P-value)",size = "Count",color = "P-value") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()+
  scale_color_manual(values = c("p < 0.001" = "#DC0000B2", "0.001 <= p < 0.01" = "#F39B7FB2", "0.01 <= p < 0.05" = "#4DBBD5B2"))



# 免疫浸润分析 ------------------------------------------------------------------
remove(list = ls()) #一键清空
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

#  ---------------- 1. Cibersort计算免疫细胞

setwd("C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\Immune")
source('Cibersort.R')

# 设置分析依赖的基础表达文件，每类免疫细胞的标志性基因及其表达，基因名字为Gene symbol
LM22.file <- "LM22.txt"
TCGA_exp.file <- "DATA.txt" 
TCGA_CIBERSORT_Results <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 1000, QN = F)
#Permutations for significance analysis是用来计算单个样本估算免疫浸润的p值，大多数文章会采用1000次。数值越大，运行时间越久

#  ---------------- 导出结果
write.csv(TCGA_CIBERSORT_Results, "TCGA_CIBERSORT_Results.csv")
write.table(TCGA_CIBERSORT_Results, "TCGA_CIBERSORT_Results.txt", sep="\t", quote =FALSE)

# -----------------比较ALN组和非ALN组

## 3. 绘图
# 3.1 数据粗处理

ALN_info <- read_csv("C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\RNA_seq\\ALN_info.csv")
colnames(ALN_info)<-c('sample','group')

library(stringr)
info <- rownames(TCGA_CIBERSORT_Results)[as.integer(substr(rownames(TCGA_CIBERSORT_Results),14,15)) == 01]
TCGA_CIBERSORT_Results <- TCGA_CIBERSORT_Results[info,]
TCGA_CIBERSORT_Results <- as.data.frame(TCGA_CIBERSORT_Results)
rownames(TCGA_CIBERSORT_Results) <- substr(rownames(TCGA_CIBERSORT_Results), 1, 12)
rownames(TCGA_CIBERSORT_Results) <- gsub("\\.","-",rownames(TCGA_CIBERSORT_Results))

TME_data <- TCGA_CIBERSORT_Results[ALN_info$sample,1:22]
TME_data$sample <- rownames(TME_data)
TME_data <- merge(TME_data, ALN_info, by = "sample", all = FALSE)

# 2.2 融合数据
TME_New = melt(TME_data)
## Using group, sample as id variables

colnames(TME_New)=c("Sample","Group","Celltype","Composition")  #设置行名
head(TME_New)

##      Group          Sample      Celltype Composition
## 1    Tumor TCGA.CV.6943.01 B cells naive 0.007651678
## 2    Tumor TCGA.CV.6959.01 B cells naive 0.019549031
## 3 Nontumor TCGA.CV.7438.11 B cells naive 0.025349204
## 4 Nontumor TCGA.CV.7242.11 B cells naive 0.032583659
## 5    Tumor TCGA.CV.7432.01 B cells naive 0.000000000
## 6 Nontumor TCGA.CV.6939.11 B cells naive 0.074282293


# 3.3 按免疫细胞占比中位数排序绘图（可选）
plot_order = TME_New[TME_New$Group=="ALN",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
write.csv(TME_New, "TME_New.csv")

# 3.3 出图
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) }

box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

box_TME;ggsave("immune_compare.pdf",box_TME,height=15,width=25,unit="cm")


# 甲基化分析 -------------------------------------------------------------------
#数据量较大，具体代码见methylation.R文件
meth_clean<-t(meth_clean)
ALN_info <- read_csv("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\ALN_info.csv")
colnames(ALN_info)<-c('sample','group')

library(stringr)
info <- rownames(meth_clean)[as.integer(substr(rownames(meth_clean),14,15)) == 01]
meth_clean <- meth_clean[info,]
meth_clean <- as.data.frame(meth_clean)
rownames(meth_clean) <- gsub("\\.","-",rownames(meth_clean))
table(substring(rownames(meth_clean),16))
Vial_A <- rownames(meth_clean)[as.character(substr(rownames(meth_clean),16,16)) == "A"]
meth_clean <- meth_clean[Vial_A,]

Vial_A_1 <- rownames(meth_clean)[grepl("A-1", rownames(meth_clean))]
meth_clean <- meth_clean[!rownames(meth_clean) %in% Vial_A_1, ]

rownames(meth_clean) <- substr(rownames(meth_clean), 1, 12)
meth_clean[6:8,1:3] 

meth_ALN_data <- meth_clean[ALN_info$sample,]
meth_ALN_data$sample <- rownames(meth_ALN_data)
meth_ALN_data <- merge(meth_ALN_data, ALN_info, by = "sample", all = FALSE)
save(meth_ALN_data, file = "TCGA-BRCA_meth_ALN.Rdata")


# LncRNA ------------------------------------------------------------------
setwd("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\")
library(readr)
library(dplyr)
filepath_gene <- "D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_lncRNA.csv"
gene_matrix <- read_csv(filepath_gene)
gene_matrix <- as.data.frame(gene_matrix)
row.names(gene_matrix) <- gene_matrix[,1]
gene_matrix <- gene_matrix[, -1]
gene_matrix[6:8,1:3] 

##1.提取出肿瘤组织
#TCGA Barcode的第14和15位代表了样本的类型，01是实体瘤，11是正常样本
table(substring(colnames(gene_matrix),14,15))
library(stringr)
tumor <- colnames(gene_matrix)[as.integer(substr(colnames(gene_matrix),14,15)) == 01]
normal <- colnames(gene_matrix)[as.integer(substr(colnames(gene_matrix),14,15)) == 11]

#将tumor样本和正常样本按顺序储存到一个矩阵中
tumor_sample <- gene_matrix[,tumor]
normal_sample <- gene_matrix[,normal]

write.csv(tumor_sample,"C:\\Users\\12923\\Desktop\\BRCA_ALN_MR\\RNA_seq\\TCGA_BRCA_tumor.csv",quote = F,row.names = T)
#exprSet_by_group <- cbind(tumor_sample,normal_sample)
#group_list <- c(rep('tumor',ncol(tumor_sample)),rep('normal',ncol(normal_sample)))

ALN_info <- read_csv("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\ALN_info.csv")
colnames(ALN_info)<-c('sample','group')

library(dplyr)
colnames(tumor_sample) <- substring(colnames(tumor_sample),1,12)
ALN_sample <- as.data.frame(ALN_info$sample[ALN_info$group=="ALN"])
Not_ALN_sample <- as.data.frame(ALN_info$sample[ALN_info$group=="Not_ALN"])

ALN_ID <- intersect(colnames(tumor_sample), ALN_sample[,1])
Not_ALN_ID <- intersect(colnames(tumor_sample), Not_ALN_sample[,1])

ALN_tumor <- t(tumor_sample[, ALN_ID])
Not_ALN_tumor <- t(tumor_sample[, Not_ALN_ID])

write.csv(ALN_tumor,"LncRNA\\ALN_LncRNA.csv",quote = F,row.names = T)
write.csv(Not_ALN_tumor,"LncRNA\\Not_ALN_LncRNA.csv",quote = F,row.names = T)

#将ALN样本和Not_ALN样本按顺序储存到一个矩阵中
exprSet_by_group <- cbind(ALN_tumor,Not_ALN_tumor)
group_list <- c(rep('ALN',ncol(ALN_tumor)),rep('Not_ALN',ncol(Not_ALN_tumor)))

#后续可直接用上述差异分析的代码


# MR数据处理 ------------------------------------------------------------------
library(readr)
library(dplyr)
filepath_mr <- "D:\\郭思蕴\\生信\\BRCA_ALN_MR\\MR_data\\final_features_MRI_91cases.csv"
MR_matrix <- read_csv(filepath_mr)
MR_matrix <- as.data.frame(MR_matrix)
row.names(MR_matrix) <- MR_matrix[,1]
MR_matrix <- MR_matrix[, -1]
MR_matrix[6:8,1:3] 
rownames(MR_matrix) <- substring(rownames(MR_matrix),1,12)

'ALN_info <- read_csv("D:\\郭思蕴\\生信\\BRCA_ALN_MR\\RNA_seq\\ALN_info.csv")
colnames(ALN_info)<-c("sample","group")

MR_ALN_matrix <- MR_matrix[intersect(rownames(MR_matrix), ALN_info$sample),]
MR_ALN_matrix$sample <- rownames(MR_ALN_matrix)
MR_ALN_matrix <- merge(MR_ALN_matrix, ALN_info, by = "sample", all = FALSE)
save(MR_ALN_matrix, file = "TCGA-BRCA_meth_ALN.Rdata")'



clinical_matrix <- as.data.frame(read_csv('D:\\郭思蕴\\生信\\BRCA_ALN_MR\\Clinical_data\\clinical_matrix.csv'))
clinical_ALN_matrix <- clinical_matrix[c('case_submitter_id', 'ajcc_pathologic_n')]
colnames(clinical_ALN_matrix)<-c('sample','ajcc_pathologic_n')

MR_ALN_matrix <- MR_matrix[intersect(rownames(MR_matrix), clinical_ALN_matrix$sample),]
MR_ALN_matrix$sample <- rownames(MR_ALN_matrix)
MR_ALN_matrix <- merge(MR_ALN_matrix, clinical_ALN_matrix, by ='sample', all = FALSE)
MR_ALN_matrix <- distinct(MR_ALN_matrix, sample,ajcc_pathologic_n, .keep_all= TRUE)

table(MR_ALN_matrix$ajcc_pathologic_n)

library(dplyr)     #数据分析包，类似于python中的pandas

ALN_i <- c('N1a','N1c', 'N2a', 'N3', 'N3a', 'N3b')
Not_ALN_i <- c('N0', 'N0 (i-)', 'N0 (i+)', 'N0(mol+)', 'N1b', 'N2b')

MR_ALN <- MR_ALN_matrix[MR_ALN_matrix$ajcc_pathologic_n %in% ALN_i, ]
MR_Not_ALN <- MR_ALN_matrix[MR_ALN_matrix$ajcc_pathologic_n %in% Not_ALN_i, ]
MR_ALN$group<-'ALN'
MR_Not_ALN$group<-'Not_ALN'

MR_ALN_data <- rbind(MR_ALN,MR_Not_ALN)
save(MR_ALN_data, file = "D:\\郭思蕴\\生信\\BRCA_ALN_MR\\MR_data\\MR_ALN_data.Rdata")
write.csv(MR_ALN_data,"D:\\郭思蕴\\生信\\BRCA_ALN_MR\\MR_data\\MR_ALN_data.csv",quote = F,row.names = F)

