setwd("/storage2/Project/CSC/RNA/03_Deconvolution/plot")

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

toolname = 'DGIST_scRNA'
mtx <- read.table(file="DGIST_result_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)

mtx <- as.matrix(mtx)
mtx_melt <- melt(mtx[,3:8])
colnames(mtx_melt) <- c('cellTypes','sample','value')

ggheatmap = ggplot(mtx_melt, aes(sample, cellTypes, fill = value)) + 
labs(title = toolname)+
geom_tile(color = "black") + 
geom_text(aes(label=round(value,2)), color='black') +
scale_fill_gradient2(low = "white", high = "red", limit = c(0, 1), space = "Lab", name = "cell type freq") + 
theme_minimal(base_size = 12, base_family = "") + 
theme(axis.text.x=element_text(angle=45,vjust=1,size=15,hjust=1),axis.text.y=element_text(vjust=1,size=15,hjust=1))+

png(paste0(toolname,"_CRITERIA_01_heatmap.png"),width=600, height=600)
ggheatmap
dev.off()
