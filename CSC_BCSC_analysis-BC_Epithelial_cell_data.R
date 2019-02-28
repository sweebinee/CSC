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
bulk = read.table(file="/storage2/Project/CSC/RNA/03_Deconvolution/CSC_RNA_TPM_rm0_HUGO_rmdup.txt", header=TRUE, stringsAsFactors=FALSE)
gene <- bulk$external_gene_name

#CIBERSORT_reference_sample_file
HP_raw_non0 <- as.matrix(HP_raw[Matrix::rowSums(HP_raw)>0,])
#write.table(HP_raw_non0,'/storage2/Project/CSC/10X/DGIST_data02/ref_sample.txt',sep = "\t", row.names=TRUE, col.names=TRUE)
> dim(HP_raw[Matrix::rowSums(HP_raw)>0,])
[1] 35126 60749
> dim(HP_raw)
[1] 36826 60749
HP_raw_non0_reduced <- HP_raw_non0[which(rownames(HP_raw_non0)%in%gene),which(colnames(HP_raw_non0) %in% cell)]
> dim(HP_raw_non0_reduced)
[1] 26008  4109
write.table(HP_raw_non0_reduced,'/storage2/Project/CSC/10X/DGIST_data02/ref_sample200_reduced.txt',sep = "\t", row.names=TRUE, col.names=TRUE)


#CIBERSORT_phenotye_file
HP_ph <- matrix(,nrow=14,ncol=60749)
rownames(HP_ph) <- c("T_cell","Monocyte","Macrophage","B_cell","Dendritic","NK_cell","GMP","CMP","Neutrophil","Haematopoietic_stem_cells","Bone_marrow_cells","Erythroblast","Myelocyte","Pro-Myelocyte")
colnames(HP_ph) <- colnames(HP_raw_non0)
for(i in 1:ncol(HP_ph)){
	cell <- colnames(HP_ph)[i]
	celltype <- HP_meta[cell,"cell_label_details"]
	HP_ph[celltype,cell] = 1
}
HP_ph[is.na(HP_ph)] <- 2
#write.table(HP_ph,'/storage2/Project/CSC/10X/DGIST_data02/pheno_sample.txt',sep = "\t", row.names=TRUE, col.names=FALSE)
#reduce tcell
tcell_list<-sample(colnames(HP_ph[,HP_ph["T_cell",]==1]),200)
B_list<-sample(colnames(HP_ph[,HP_ph["B_cell",]==1]),200)
DC_list<-sample(colnames(HP_ph[,HP_ph["Dendritic",]==1]),200)                   
Macro_list<-sample(colnames(HP_ph[,HP_ph["Macrophage",]==1]),200)                   
Mono_list<-sample(colnames(HP_ph[,HP_ph["Monocyte",]==1]),200)                   
Neutro_list<-sample(colnames(HP_ph[,HP_ph["Neutrophil",]==1]),200)                   
NK_list<-sample(colnames(HP_ph[,HP_ph["NK_cell",]==1]),200)                   

BM_list<-colnames(HP_ph[,HP_ph["Bone_marrow_cells",]==1])
CMP_list<-colnames(HP_ph[,HP_ph["CMP",]==1])
Erythro_list<-colnames(HP_ph[,HP_ph["Erythroblast",]==1])
GMP_list<-colnames(HP_ph[,HP_ph["GMP",]==1])
Haema_list<-colnames(HP_ph[,HP_ph["Haematopoietic_stem_cells",]==1])
Myelo_list<-colnames(HP_ph[,HP_ph["Myelocyte",]==1])
ProM_list<-colnames(HP_ph[,HP_ph["Pro-Myelocyte",]==1])

cell <- union(union(union(union(union(union(union(union(union(union(union(union(union(tcell_list,B_list),DC_list),Macro_list),Mono_list),Neutro_list),NK_list),BM_list),CMP_list),Erythro_list),GMP_list),Haema_list),Myelo_list),ProM_list)
#4109 cells : 500
#2709 cells : 300
#1309 cells : 100
#2009 cells : 200

HP_ph_reduced<-HP_ph[,which(colnames(HP_ph) %in% cell)]
write.table(HP_ph_reduced,'/storage2/Project/CSC/10X/DGIST_data02/pheno_sample200_reduced.txt',sep = "\t", row.names=TRUE, col.names=FALSE)

for(i in 1:10){
	tcell_list<-sample(colnames(HP_ph[,HP_ph["T_cell",]==1]),200)
	B_list<-sample(colnames(HP_ph[,HP_ph["B_cell",]==1]),200)
	DC_list<-sample(colnames(HP_ph[,HP_ph["Dendritic",]==1]),200)                   
	Macro_list<-sample(colnames(HP_ph[,HP_ph["Macrophage",]==1]),200)                   
	Mono_list<-sample(colnames(HP_ph[,HP_ph["Monocyte",]==1]),200)                   
	Neutro_list<-sample(colnames(HP_ph[,HP_ph["Neutrophil",]==1]),200)                   
	NK_list<-sample(colnames(HP_ph[,HP_ph["NK_cell",]==1]),200)                   
	BM_list<-colnames(HP_ph[,HP_ph["Bone_marrow_cells",]==1])
	CMP_list<-colnames(HP_ph[,HP_ph["CMP",]==1])
	Erythro_list<-colnames(HP_ph[,HP_ph["Erythroblast",]==1])
	GMP_list<-colnames(HP_ph[,HP_ph["GMP",]==1])
	Haema_list<-colnames(HP_ph[,HP_ph["Haematopoietic_stem_cells",]==1])
	Myelo_list<-colnames(HP_ph[,HP_ph["Myelocyte",]==1])
	ProM_list<-colnames(HP_ph[,HP_ph["Pro-Myelocyte",]==1])
	cell <- union(union(union(union(union(union(union(union(union(union(union(union(union(tcell_list,B_list),DC_list),Macro_list),Mono_list),Neutro_list),NK_list),BM_list),CMP_list),Erythro_list),GMP_list),Haema_list),Myelo_list),ProM_list)
	HP_raw_non0_reduced <- HP_raw_non0[which(rownames(HP_raw_non0)%in%gene),which(colnames(HP_raw_non0) %in% cell)]
	write.table(HP_raw_non0_reduced,paste0('/storage2/Project/CSC/10X/DGIST_data02/CIBERSORT_200/ref_sample200_reduced_',i,'.txt'),sep = "\t", row.names=TRUE, col.names=TRUE)
	HP_ph_reduced<-HP_ph[,which(colnames(HP_ph) %in% cell)]
	write.table(HP_ph_reduced,paste0('/storage2/Project/CSC/10X/DGIST_data02/CIBERSORT_200/pheno_sample200_reduced_',i,'.txt'),sep = "\t", row.names=TRUE, col.names=FALSE)
}


