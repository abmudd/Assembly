#!/usr/bin/env Rscript

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }

args = strsplit(commandArgs(TRUE),' ');
data <- read.csv(args[[1]], header=FALSE)
cnames <- data[1,2:ncol(data)]
rnames <- data[2:nrow(data),1]
mat_data <- data.matrix(data[2:nrow(data),2:ncol(data)])
rownames(mat_data) <- rnames
colnames(mat_data) <- cnames

my_palette <- colorRampPalette(c("black","black","blue","white","red"))(n = 299)
png(paste(args[[1]],"png",sep='.'), width = 10*300, height = 10*300, res = 300, pointsize = 15)
heatmap.2(mat_data, Rowv=NA, Colv=NA, main=" ", col=my_palette, xlab="Fissions", ylab="Fusions", scale="none", margins=c(5,5), cexCol=0.75, cexRow=0.75, trace="none", density.info="none", symkey=FALSE, dendrogram="none", key.par=list(mar=c(5,1,5,0)), lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2, 6), lwid=c(1, 10, 1))
#title(ylab=paste("Fusions"), line = 0)
#title(main=paste(args[[1]],"Heat Map"), line = 2.5)
dev.off()
