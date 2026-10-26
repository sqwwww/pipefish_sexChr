library(qqman)
library(ggplot2)

args=commandArgs(T)
fileName=args[1]
outName_prefix=args[2]
mydata <- read.table(file = fileName, header = F, sep = " ")
#mydata <- read.table(file = "test.txt", header = F, sep = " ")
#head(mydata)
names(mydata) <- c("SNP","CHR","BP","P")
#head(mydata)
source("myManhattanFunction.R")
p <- myManhattan(mydata, graph.title = "My Manhattan Plot",suggestiveline = FALSE,genomewideline = 1e-5,highlight = 1e-5)
ggsave(filename = paste0("MH_",outName_prefix,".jpeg"), plot = p, device = "jpeg")
jpeg(file=paste0("qqplot_",outName_prefix,".jpeg"))
qq(mydata$P)
dev.off()
