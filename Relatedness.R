library(ggplot2)
pdf("relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
ggplot(data=relatedness,aes(x=Z0,y=Z1,color=RT,alpha=0.5))+geom_point()

pdf("zoom_relatedness.pdf")
relatedness_zoom = read.table("zoom_pihat.genome", header=T)
ggplot(data=relatedness_zoom,aes(x=Z0,y=Z1,color=RT,alpha=0.5))+geom_point()

pdf("hist_relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
dev.off() 

