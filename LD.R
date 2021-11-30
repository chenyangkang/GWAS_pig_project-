library(LDheatmap)
library(genetics)
pdf("LD.pdf")
for(i in 1:18){
  print(i)
  name <- read.table(paste0('d1_6_chr',i,'.map'))
  mx <- as.matrix(read.table(paste0('d1_6_chr',i,'.ld')))
  colnames(mx) <- name[,2]
  rownames(mx) <- name[,2]
  dist<-name[,4]
  #作图
  ld <- LDheatmap(mx,genetic.distances = dist)
  color = c(colorRampPalette(colors = c("red","white"))(80),colorRampPalette(colors = c("white","grey"))(20))
  color2 <- colorRampPalette(c("red","#FA9FB5","#BDBDBD"))
  LDheatmap(ld,color = color2(100),SNP.name = c('SNP1','SNP2'))
}
dev.off()

