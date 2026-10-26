#library(qqman)
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
p <- myManhattan(mydata, graph.title = "My Manhattan Plot",suggestiveline = 4.09e-7,genomewideline = 2.05e-8,highlight = 2.05e-8,highlight.col = "#e25822",col = c("#a1caf1", "#0067a5"))
ggsave(filename = paste0("MH_",outName_prefix,".png"), plot = p, device = "png", width = 5, height = 8)
#jpeg(file=paste0("qqplot_",outName_prefix,".jpeg"))
#qq(mydata$P)
#dev.off()
