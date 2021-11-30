args=commandArgs(T)
#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman")  
#results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
#jpeg("Logistic_manhattan.jpeg")
#manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")
#dev.off()

results_as <- read.table(args[1], head=TRUE)
pdf(paste0(args[2],".manhatton.pdf"))
manhattan(results_as,chr="CHR",bp="POS",p="P",snp="SNP", main = paste0("Manhattan plot:",args[2]))
dev.off()




