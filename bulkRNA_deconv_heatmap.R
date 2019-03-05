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

###############################################################################
#scatter plot
library(plotly)

mtx <- read.table(file="DGIST_result_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)
mtx <- as.matrix(mtx)
mtx_melt <- melt(mtx[,3:8])
#> mtx_melt
#                        Var1 Var2        value
#1                     B_cell PE24 0.0610323672
colnames(mtx_melt)[3] <- "scRNA"

toolnames = c('MCPcounter','CIBERSORT','Xcell')
for(i in toolnames){
  assign(i, read.delim(paste0("01CRITERIA_",i,".txt"), header=T, stringsAsFactors=F, row.names=1))
  assign(paste0(i,"_melt"),melt(as.matrix(get(i))))
  assign(paste0(i,"_merge"),merge(mtx_melt,get(paste0(i,"_melt")),by=c('Var1','Var2')))
  png(paste0(i,"_CRITERIA_01_dot.png"),width=600, height=600)
  plot(x=get(paste0(i,"_merge"))$scRNA,y=get(paste0(i,"_merge"))$value)
  dev.off()
}

colnames(MCPcounter_merge)[4] <- "MCPcounter"
merge_mtx <- merge(MCPcounter_merge,CIBERSORT_melt,by=c('Var1','Var2'))
colnames(merge_mtx)[5]<-"CIBERSORT"
merge_mtx <- merge(merge_mtx,Xcell_melt,by=c('Var1','Var2'))
colnames(merge_mtx)[6]<-"Xcell"

# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
png("CRITERIA_01_dot.png",width=600, height=600)
pairs(merge_mtx[,3:6], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="Title")
dev.off()


#scRNA-CIBERSORT, scRNA-Xcell only scatter plot
library(ggrepel)

for(i in 1:nrow(merge_mtx)){
	rownames(merge_mtx)[i]<-paste0(merge_mtx$Var1[i],":",merge_mtx$Var2[i])
}

merge_mtx<-merge_mtx[-c(79:84),]

# Add text to the plot
.labs <- rownames(merge_mtx)
b <- ggplot(merge_mtx, aes(x = scRNA, y = CIBERSORT))

png("CRITERIA_03_dot_scRNA_CIBERSORT_noT.png",width=800, height=600)
b + xlim(0.001, 1)+ylim(0.001, 1)+
  geom_point(aes(color = Var1)) +
  geom_smooth(method='lm',se = FALSE, fullrange = TRUE)+
  ggpubr::stat_cor(label.x = 0.003)+
  geom_text_repel(aes(label = .labs,  color = Var1), size = 3)+
  scale_color_manual(values = c("#FF0000", "#FF5500","#FFAA00","#FFFF00","#AAFF00","#00FF2B","#00FFD4","#00D4FF","#00AAFF","#0055FF","#5500FF","#AA00FF","#FF00AA","#FF0055"))+
  labs(title = "scRNA vs CIBERSORT\n", color = "cellTypes\n")
dev.off()



