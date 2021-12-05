setwd("C:/Users/min/Desktop/统计生物与网络生物学/富集分析")
getwd()

library(rtracklayer)

#####################################################################################################################################
expr_df_T3 <- read.csv("T3_result_df_gene.csv")
head(expr_df_T3)

###基因ID转换
library(org.Ss.eg.db)
library(AnnotationDbi)
help("org.Ss.eg.db")
columns(org.Ss.eg.db)
keys(org.Ss.eg.db, keytype="ACCNUM") %>% head()

expr_df_T3$entrez <- mapIds(org.Ss.eg.db,
                     keys=expr_df_T3$gene_id,
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")
write.csv(expr_df_T3,"T3_result_df_gene_ENTREZIDID.csv")

library(clusterProfiler)
####GO分析
GO_database <- 'org.Ss.eg.db'
GO<-enrichGO( expr_df_T3$entrez,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.5,#设定p值阈值
              qvalueCutoff = 0.5,#设定q值阈值
              readable = T)
write.csv(GO,"T4_result_df_gene_GO分析_MF-BP-CC.csv")
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
##富集基因与所在功能集/通路集的关联网络图：
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::heatplot(GO,showCategory = 50)#基因-通路关联热图

###KEGG分析
library(enrichplot)
library(GOSemSim)
library(DOSE)
kk <- enrichKEGG(gene = expr_df_T3$entrez,
                 organism ="ssc",
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 0.5,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)
write.csv(kk,"T3_result_df_gene_KEGG分析全部结果.csv")
head(kk)[,1:6]
###作图
barplot(kk)
dotplot(kk)

###富集基因与所在功能集/通路集的关联网络图：
enrichplot::cnetplot(kk,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
kk2 <- pairwise_termsim(kk)
enrichplot::emapplot(kk2, color = "p.adjust", layout = "kk")\

##########################################################################################################################
####韦恩图：
###VennDiagram函数包最大能绘制5个数据集合的韦恩图
library(VennDiagram)
library(futile.logger)
dat <- read.table('T3_result_df_genes_韦恩图数据准备【算法汇总】.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(dat)
venn_list <- list(Closeness = dat$Closeness, Degree = dat$Degree,DMNC=dat$DMNC,MCC=dat$MCC)
venn.diagram(venn_list, filename = 'T4_result_df_genes_DEGs-cytoHubba selection.png', imagetype = 'png', 
             fill = c('red', 'blue', 'green', 'orange'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green', 'orange'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green', 'orange'), cex = 1.5, fontfamily = 'serif')
####函数参数解释：
#cex = 2,        #区域内部数字的字体大小 
#cat.cex = 2,    #分类名称的字体大小 
#cat.fontfamily = "serif",     #分类字体
#cat.col = c("darkblue", "darkgreen", "orange", "grey50"),  #分类颜色
####获得各组之间的交集元素
#继续以上述4个分组为例，组间交集元素获得
inter <- get.venn.partitions(venn_list)
head(inter)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(5, 6)], 'T3_result_df_genes_算法筛选的Hub基因交集-韦恩图筛选.txt', row.names = FALSE, sep = '\t', quote = FALSE)


###########################T4分析
expr_df_T4 <- read.csv("T4_result_df_gene.csv")
head(expr_df_T4)

###基因ID转换
library(org.Ss.eg.db)
library(AnnotationDbi)
help("org.Ss.eg.db")
columns(org.Ss.eg.db)
keys(org.Ss.eg.db, keytype="ACCNUM") %>% head()

expr_df_T4$entrez <- mapIds(org.Ss.eg.db,
                            keys=expr_df_T4$gene_id,
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
write.csv(expr_df_T4,"T4_result_df_gene_ENTREZIDID.csv")

library(clusterProfiler)
####GO分析
GO_database <- 'org.Ss.eg.db'
GO1<-enrichGO( expr_df_T4$entrez,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.5,#设定p值阈值
              qvalueCutoff = 0.5,#设定q值阈值
              readable = T)
write.csv(GO1,"T4_result_df_gene_GO分析_MF-BP-CC.csv")
barplot(GO1, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
dotplot(GO1, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
##富集基因与所在功能集/通路集的关联网络图：
enrichplot::cnetplot(GO1,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::heatplot(GO1,showCategory = 50)#基因-通路关联热图

###KEGG分析
library(enrichplot)
library(GOSemSim)
library(DOSE)
kk1 <- enrichKEGG(gene = expr_df_T4$entrez,
                 organism ="ssc",
                 pvalueCutoff = 0.5,
                 qvalueCutoff = 0.5,
                 minGSSize = 1,
                 #readable = TRUE ,
                 use_internal_data =FALSE)
write.csv(kk1,"T4_result_df_gene_KEGG分析全部结果.csv")
head(kk1)[,1:6]
###作图
barplot(kk1)
dotplot(kk1)

###富集基因与所在功能集/通路集的关联网络图：
enrichplot::cnetplot(kk1,circular=FALSE,colorEdge = TRUE)#circluar为指定是否环化，基因过多时建议设置为FALSE
kk3 <- pairwise_termsim(kk1)
enrichplot::emapplot(kk3, color = "p.adjust", layout = "kk")\

##########################################################################################################################
####韦恩图：
###VennDiagram函数包最大能绘制5个数据集合的韦恩图
library(VennDiagram)
library(futile.logger)
dat1 <- read.table('T4_result_df_genes_韦恩图数据准备【算法汇总】.txt', header = TRUE, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
head(dat1)
venn_list1 <- list(Closeness = dat$Closeness, Degree = dat$Degree,DMNC=dat$DMNC,MCC=dat$MCC)
venn.diagram(venn_list1, filename = 'T4_result_df_genes_DEGs-cytoHubba selection.png', imagetype = 'png', 
             fill = c('red', 'blue', 'green', 'orange'), alpha = 0.50, 
             cat.col = c('red', 'blue', 'green', 'orange'), cat.cex = 1.5, cat.fontfamily = 'serif',
             col = c('red', 'blue', 'green', 'orange'), cex = 1.5, fontfamily = 'serif')
####函数参数解释：
#cex = 2,        #区域内部数字的字体大小 
#cat.cex = 2,    #分类名称的字体大小 
#cat.fontfamily = "serif",     #分类字体
#cat.col = c("darkblue", "darkgreen", "orange", "grey50"),  #分类颜色
####获得各组之间的交集元素
#继续以上述4个分组为例，组间交集元素获得
inter1 <- get.venn.partitions(venn_list1)
head(inter1)
for (i in 1:nrow(inter1)) inter1[i,'values'] <- paste(inter1[[i,'..values..']], collapse = ', ')
write.table(inter1[-c(5, 6)], 'T4_result_df_genes_算法筛选的Hub基因交集-韦恩图筛选.txt', row.names = FALSE, sep = '\t', quote = FALSE)
