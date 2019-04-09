start at April 08 

library(Seurat)
library(ggplot2)
library(RColorBrewer)

setwd('/storage2/Project/CSC/10X/DGIST_data02')

#LM22 marker gene
CIBERSORT_gene <- read.table('/storage2/Project/CSC/10X/DGIST_data02/test_cellType/CIBERSORT_DEG_list.txt', header=TRUE, row.names=1, sep='\t')
for(i in colnames(CIBERSORT_gene)){
	assign(paste0(i,'_DEG'),rownames(CIBERSORT_gene[CIBERSORT_gene[,i]==1,]))
}
#cluster marker gene
clusters <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")
for(i in clusters){
  assign(paste0("DEG_",i), readRDS(paste0("DEG_",i,".rds")))
}
for(i in clusters){
  deg <- get(paste0("DEG_",i))
  assign(paste0("DEG_cut_",i), deg[deg$p_val_adj<0.05&abs(deg$avg_logFC)>1,])
  assign(paste0("DEG_cut_list_",i),rownames(get(paste0("DEG_cut_",i))))
}

#intersect table(row=LM22, col=Clusters)
Int_table <- matrix(,ncol=22,nrow=22)
colnames(Int_table)=clusters
rownames(Int_table)=colnames(CIBERSORT_gene)
for(i in colnames(Int_table)){
	cluster_deg <- get(paste0("DEG_cut_list_",i))
	for(j in rownames(Int_table)){
		cell_deg <- get(paste0(j,'_DEG'))
		INT <- length(intersect(cluster_deg,cell_deg))
		UNI <- length(union(cluster_deg,cell_deg))
		Int_table[j,i] <- round(INT/UNI,3)*100
	}
}

mtx_melt <- melt(Int_table)
colnames(mtx_melt)<-c("cellTypes","Clusters","value")

ggheatmap = ggplot(mtx_melt, aes(Clusters, cellTypes, fill = value)) + 
labs(title = "weighted CellType_DEG & Cluster_DEG \n intersect / union")+
geom_tile(color = "black") + 
geom_text(aes(label=round(value,2)), color='black') +
scale_fill_gradient2(low = "white", high = "red", limit = c(0, 12), space = "Lab", name = ".") + 
theme_minimal(base_size = 12, base_family = "") + 
theme(axis.text.x=element_text(angle=45,vjust=1,size=15,hjust=1),axis.text.y=element_text(vjust=1,size=15,hjust=1))

png("intersect_table_heatmap_weighted.png",width=800, height=600)
ggheatmap
dev.off()
################################################################
##################################################################
#with weighting
count_CIBERSORT_DEG <- ''
for(i in rownames(CIBERSORT_gene)){
	count_CIBERSORT_DEG[i] <- as.character(rowSums(CIBERSORT_gene[i,]))
}
> table(count_CIBERSORT_DEG)
count_CIBERSORT_DEG
      1   2   3   4   5   6   7   8   9 
    262 146  58  42  13  11  10   2   3 

#intersect table(row=LM22, col=Clusters)
Int_table <- matrix(,ncol=22,nrow=22)
colnames(Int_table)=clusters
rownames(Int_table)=colnames(CIBERSORT_gene)
for(i in colnames(Int_table)){
	cluster_deg <- get(paste0("DEG_cut_list_",i))
	for(j in rownames(Int_table)){
		cell_deg <- get(paste0(j,'_DEG'))
		INT_score <- intersect(cluster_deg,cell_deg)
		iscore = 0
		for(x in INT_score){
			iscore = iscore + 1/int(count_CIBERSORT_DEG[x])
			print(iscore)
		}
		UNI_score <- union(cluster_deg,cell_deg)
		Int_table[j,i] <- round(INT/UNI,3)*100
	}
}
################################################################
##################################################################
cellType = "Dendritic.cells.resting" 
cluster = 17
CIBERSORT <-get(paste0(cellType,"_DEG"))
CSC <- get(paste0("DEG_cut_list_",cluster))

INT <- intersect(CIBERSORT,CSC)

length(CIBERSORT)
length(CSC)
length(INT)

diff <- setdiff(CSC,INT)

for(target in diff){
#pdf(paste0("Neutrophil_",target,".pdf"))
cols <- c('#D5D8DC',brewer.pal(9, "Reds"))
df = data.frame(x=blood@reductions$tsne@cell.embeddings[, "tSNE_1"], 
                y=blood@reductions$tsne@cell.embeddings[, "tSNE_2"], 
                expression=blood@assays$RNA@scale.data[target,])
g<-ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=1) + 
  scale_colour_gradientn(colours = cols ) +
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=20), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  )+ ggtitle(paste0("B.cells.memory's new marker : ",target))
ggsave(paste0("B.cells.memory_",target,".png"),g)
#dev.off()
}