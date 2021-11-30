args=commandArgs(T)
#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
#results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
#pdf("T3-QQ-Plot_logistic.pdf")
#qq(results_log$P, main = "T3 Q-Q plot of GWAS p-values : log")
#dev.off()

results_as <- read.table(args[1], head=TRUE)
pdf(paste0(args[2],".qqplot.pdf"))
qq(results_as$P, main = paste0(args[2],"-Q-Q plot of GWAS p-values : log"))
dev.off()

