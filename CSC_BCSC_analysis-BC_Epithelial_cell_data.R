start on 25Feb
BC_Epithelial_cell_data

library(Seurat)
library(ggplot2)
library(RColorBrewer)

setwd('/storage2/Project/CSC/10X/DGIST_data02')
DGIST=readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Epithelial_cell_data/EpithelialCells_SeuratSet.rds")
#BCSC_meta$ = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Epithelial_cell_data/Epithelial_metadata.rds")
#BCSC_raw = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Epithelial_cell_data/Epithelial_Rawcount.rds")
#BCSC_norm = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Epithelial_cell_data/Epithelial_normalizedExprs.rds")
DGIST@assays        DGIST@active.ident  DGIST@reductions    DGIST@version       
DGIST@meta.data     DGIST@graphs        DGIST@project.name  DGIST@commands      
DGIST@active.assay  DGIST@neighbors     DGIST@misc          DGIST@tools 

pdf("DGIST_data02_BCSC_tsne.pdf")
DimPlot(DGIST,reduction = "tsne")
dev.off()

#윤정교 교수님 wnt signal
target = 'RHF43'

pdf(paste0("dr.Yoon_wnt_",target,".pdf"))
cols <- brewer.pal(9,"YlOrRd")
df = data.frame(x=DGIST@reductions$tsne@cell.embeddings[, "tSNE_1"], 
                y=DGIST@reductions$tsne@cell.embeddings[, "tSNE_2"], 
                expression=DGIST@assays$RNA@scale.data[target,])
ggplot(df,aes(x=x, y=y, colour=expression)) + 
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
  )+ ggtitle(paste0("Wnt signal : ",target))
dev.off()

#bulkRNA-seq deconvolution 비교
#cluster별 variable gene 찾고 CIBERSORT, Xcell에 대입
HP_raw = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_Rawcount.rds")
HP_norm = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_normalizedExprs.rds")
HP_meta = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/Hematopoietic_metadata.rds")
table(HP_meta$cell_label_details)
                   B_cell         Bone_marrow_cells                       CMP 
                     3017                        17                        10 
                Dendritic              Erythroblast                       GMP 
                     4322                         6                        34 
Haematopoietic_stem_cells                Macrophage                  Monocyte 
                      191                      4673                      8461 
                Myelocyte                Neutrophil                   NK_cell 
                      339                      1183                      3905 
            Pro-Myelocyte                    T_cell 
                       12                     34579 

#다쓰면 데이터 너무커서 CIBERSORT에 안올라감 1/8 = 7593 CELL 정도 써야함 (500개이상이면 많은거) 
#300개 random sampling
#bulkRNA에서 쓴 gene만 보자
#bulk = read.table(file="/storage2/Project/CSC/RNA/03_Deconvolution/CSC_RNA_TPM_rm0_HUGO_rmdup.txt", header=TRUE, stringsAsFactors=FALSE)
#gene <- bulk$external_gene_name

#CIBERSORT_reference_sample_file
HP_norm_non0 <- as.matrix(HP_norm[Matrix::rowSums(HP_norm)>0,])
#HP_norm_non0_reduced <- HP_norm_non0[which(rownames(HP_norm_non0)%in%gene),]

CIBERSORT_gene <- read.table('/storage2/Project/CSC/10X/DGIST_data02/CIBERSORT_signature_gene_list.txt')$V1
S12_gene <- read.table('/storage2/Project/CSC/10X/DGIST_data02/S12_signature_gene_list.txt')$V1
S3_gene <- read.table('/storage2/Project/CSC/10X/DGIST_data02/S3_signature_gene_list.txt')$V1
sig_gene <- union(union(CIBERSORT_gene,S12_gene),S3_gene)

HP_norm_non0_reduced <- HP_norm_non0[which(rownames(HP_norm_non0)%in%union_DEG),]

HP_celltype_exp <- matrix(, nrow=nrow(HP_norm_non0_reduced), ncol=112)
sample=unique(HP_meta$sample_name)
celltypes <- c("T_cell","Monocyte","Macrophage","B_cell","Dendritic","NK_cell","GMP","CMP","Neutrophil","Haematopoietic_stem_cells","Bone_marrow_cells","Erythroblast","Myelocyte","Pro-Myelocyte")
coln <- c()
for(i in celltypes){
	for(j in sample){
		x <- paste0(i,":",j)
		coln <- append(coln,x)
	}
}

colnames(HP_celltype_exp) <- coln 
rownames(HP_celltype_exp) <- rownames(HP_norm_non0_reduced)

for(i in celltypes){
	for(j in sample){
		cell <- rownames(HP_meta[HP_meta$cell_label_details==i&HP_meta$sample_name==j,])
		for(z in rownames(HP_celltype_exp)){
			avg <- mean(HP_norm_non0_reduced[z,which(colnames(HP_norm_non0_reduced)%in%cell)],na.rm = TRUE)
			HP_celltype_exp[z,paste0(i,":",j)] <- avg
			print(paste0(i,":",j))
		}
	}
}

dim(HP_celltype_exp)
HP_celltype_exp[is.na(HP_celltype_exp)]<-0

write.table(HP_celltype_exp,'/storage2/Project/CSC/10X/DGIST_data02/ref_sample_reduced.txt',sep = "\t", row.names=TRUE, col.names=TRUE)

#weighting
HP_celltype_exp <- read.table("/storage2/Project/CSC/10X/DGIST_data02/mark05_ref.txt",sep='\t',header=TRUE,row.names=1)


#CIBERSORT_phenotye_file
HP_ph <- matrix(,nrow=14,ncol=112)
rownames(HP_ph) <- c("T_cell","Monocyte","Macrophage","B_cell","Dendritic","NK_cell","GMP","CMP","Neutrophil","Haematopoietic_stem_cells","Bone_marrow_cells","Erythroblast","Myelocyte","Pro-Myelocyte")
colnames(HP_ph) <- colnames(HP_celltype_exp)

for(i in colnames(HP_ph))){
	celltype <- strsplit(i,':')[[1]][1]
	HP_ph[celltype,i] = 1
}
HP_ph[is.na(HP_ph)] <- 2
write.table(HP_ph,'/storage2/Project/CSC/10X/DGIST_data02/pheno_sample.txt',sep = "\t", row.names=TRUE, col.names=FALSE)


#scRNA-seq DEG
blood = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/BC_blood_seurat.rds")
blood_cellTypes = blood
ident_mtx <- as.matrix(blood_cellTypes@active.ident)
for(i in rownames(ident_mtx)){
	ident_mtx[i,1]<-HP_meta[i,"cell_label_details"]
}
blood_cellTypes@active.ident<-factor(ident_mtx[,1])
> levels(blood_cellTypes@active.ident)
 [1] "B_cell"                    "Bone_marrow_cells"        
 [3] "CMP"                       "Dendritic"                
 [5] "Erythroblast"              "GMP"                      
 [7] "Haematopoietic_stem_cells" "Macrophage"               
 [9] "Monocyte"                  "Myelocyte"                
[11] "Neutrophil"                "NK_cell"                  
[13] "Pro-Myelocyte"             "T_cell"    

for(i in levels(blood_cellTypes@active.ident)){
	assign(paste0(i,"_DEG"), FindMarkers(blood_cellTypes, ident.1=i, ident.2=NULL, only.pos=TRUE))
}

for(i in levels(blood_cellTypes@active.ident)){
	saveRDS(get(paste0(i,"_DEG")),paste0(i,"_DEG.rds"))
}


celltypes <- c("T_cell","Monocyte","Macrophage","B_cell","Dendritic","NK_cell","GMP","CMP","Neutrophil","Haematopoietic_stem_cells","Bone_marrow_cells","Erythroblast","Myelocyte","ProMyelocyte")
for(i in celltypes){
  assign(paste0(i,"_DEG"), readRDS(paste0(i,"_DEG.rds")))
}

for(i in celltypes){
  deg <- get(paste0(i,"_DEG"))
  assign(paste0(i,"_DEG_cut"), deg[deg$p_val_adj<0.05,])
  degList <- 
  assign(paste0(i,"_DEG_cut_list"),rownames(get(paste0(i,"_DEG_cut"))))
}

union_DEG<-union(union(union(union(union(union(union(union(union(union(union(union(union(T_cell_DEG_cut_list,Monocyte_DEG_cut_list),Macrophage_DEG_cut_list),B_cell_DEG_cut_list),Dendritic_DEG_cut_list),NK_cell_DEG_cut_list),GMP_DEG_cut_list),CMP_DEG_cut_list),Neutrophil_DEG_cut_list),Haematopoietic_stem_cells_DEG_cut_list),Bone_marrow_cells_DEG_cut_list),Erythroblast_DEG_cut_list),Myelocyte_DEG_cut_list),ProMyelocyte_DEG_cut_list)

#make weight table
DEG_w <- matrix(,nrow=length(union_DEG), ncol=1)
rownames(DEG_w) <- union_DEG
for(i in union_DEG){
  w = 0
  for(j in celltypes){
    deg <- get(paste0(j,"_DEG_cut_list"))
    if(i %in% deg) w = w+1
  }
  DEG_w[i,1] <- w
}
> table(DEG_w)
DEG_w
   1    2    3    4    5    6    7 
1683  812  463  190   64   49    3 


