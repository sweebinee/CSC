
library(Seurat)
library(ggplot2)
library(RColorBrewer)

setwd('/storage2/Project/CSC/10X/DGIST_data02')

#scRNA-seq DEG by cluster 
blood = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/BC_blood_seurat.rds")
HP_raw = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_Rawcount.rds")
#HP_norm = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_normalizedExprs.rds")
HP_meta = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_metadata.rds")

> levels(blood@active.ident)
  [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
[16] "15" "16" "17" "18" "19" "20" "21"
clusters <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")

for(i in levels(blood@active.ident)){
	assign(paste0("DEG_",i), FindMarkers(blood, ident.1=i, ident.2=NULL,only.pos=TRUE))
}

for(i in levels(blood@active.ident)){
	saveRDS(get(paste0("DEG_",i)),paste0("DEG_onlyPOS",i,".rds"))
}

for(i in clusters){
  assign(paste0("DEG_",i), readRDS(paste0("DEG_",i,".rds")))
}

for(i in clusters){
  deg <- get(paste0("DEG_",i))
  assign(paste0("DEG_cut_",i), deg[deg$p_val_adj<0.05,])
  assign(paste0("DEG_cut_list_",i),rownames(get(paste0("DEG_cut_",i))))
}

union_DEG<-union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(union(DEG_cut_list_0,DEG_cut_list_1),DEG_cut_list_2),DEG_cut_list_3),DEG_cut_list_4),DEG_cut_list_5),DEG_cut_list_6),DEG_cut_list_7),DEG_cut_list_8),DEG_cut_list_9),DEG_cut_list_10),DEG_cut_list_11),DEG_cut_list_12),DEG_cut_list_13),DEG_cut_list_14),DEG_cut_list_15),DEG_cut_list_16),DEG_cut_list_17),DEG_cut_list_18),DEG_cut_list_19),DEG_cut_list_20),DEG_cut_list_21)
> length(union_DEG)
[1] 5017

png("heatmap_cluster.png")
DoHeatmap(blood,genes.use = union_DEG,slim.col.label = TRUE,remove.key = TRUE,cex.row=0.1)
dev.off()


#################
> table(HP_meta$clusters)

   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
7001 6940 5721 4655 4580 4518 4340 4172 3630 2672 2112 1884 1653 1246 1186 1169 
  16   17   18   19   20   21 
1077  923  480  383  296  111 

HP_raw_non0 <- as.matrix(HP_raw[Matrix::rowSums(HP_raw)>0,])
HP_raw_non0_reduced <- HP_raw_non0[which(rownames(HP_raw_non0)%in%union_DEG),]
> dim(HP_raw_non0_reduced)
[1]  5008 60749

#make ref_sample.txt
HP_celltype_exp <- matrix(, nrow=nrow(HP_raw_non0_reduced), ncol=22)
colnames(HP_celltype_exp) <- clusters 
rownames(HP_celltype_exp) <- rownames(HP_raw_non0_reduced)

for(i in clusters){
    cell <- rownames(HP_meta[HP_meta$clusters==i,])
    for(z in rownames(HP_celltype_exp)){
      avg <- mean(HP_raw_non0_reduced[z,which(colnames(HP_raw_non0_reduced)%in%cell)],na.rm = TRUE)
      HP_celltype_exp[z,i] <- avg
      print(i)
    }
}

write.table(HP_celltype_exp,'/storage2/Project/CSC/10X/DGIST_data02/ref_cluster_logFC1.txt',sep = "\t", row.names=TRUE, col.names=TRUE)


###
#DGIST cluster - PEXX table
cell_count<-matrix(,nrow=22,ncol=8)
sample=unique(HP_meta$sample_name)
colnames(cell_count)<-sample
rownames(cell_count)<-clusters
for(i in clusters){
  for(j in sample){
    cell <- length(rownames(HP_meta[HP_meta$clusters==i&HP_meta$sample_name==j,]))
    cell_count[i,j]<-cell
  }
}
write.table(cell_count,'/storage2/Project/CSC/10X/DGIST_data02/DGIST_cluster_sample_count.txt',sep = "\t", row.names=TRUE, col.names=TRUE)

