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
HP_raw_non0_reduced <- HP_raw_non0[which(rownames(HP_raw_non0)%in%gene),]

HP_celltype_exp <- matrix(, nrow=nrow(HP_raw_non0_reduced), ncol=112)
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
rownames(HP_celltype_exp) <- rownames(HP_raw_non0_reduced)

for(i in celltypes){
	for(j in sample){
		cell <- rownames(HP_meta[HP_meta$cell_label_details==i&HP_meta$sample_name==j,])
		for(z in rownames(HP_celltype_exp)){
			avg <- mean(HP_raw_non0_reduced[z,which(colnames(HP_raw_non0_reduced)%in%cell)])
			HP_celltype_exp[z,paste0(i,":",j)] <- avg
			print(paste0(i,":",j))
		}
	}
}
> dim(HP_celltype_exp)
[1] 26008    14
write.table(HP_celltype_exp,'/storage2/Project/CSC/10X/DGIST_data02/ref_sample_reduced.txt',sep = "\t", row.names=TRUE, col.names=TRUE)


#CIBERSORT_phenotye_file
HP_ph <- matrix(,nrow=14,ncol=112)
rownames(HP_ph) <- c("T_cell","Monocyte","Macrophage","B_cell","Dendritic","NK_cell","GMP","CMP","Neutrophil","Haematopoietic_stem_cells","Bone_marrow_cells","Erythroblast","Myelocyte","Pro-Myelocyte")
colnames(HP_ph) <- colnames(HP_celltype_exp)
for(i in colnames(HP_ph))){
	celltype <- strsplit(i,':')[1]
	HP_ph[celltype,i] = 1
}
colsum(HP_ph)

HP_ph[is.na(HP_ph)] <- 2
#write.table(HP_ph,'/storage2/Project/CSC/10X/DGIST_data02/pheno_sample.txt',sep = "\t", row.names=TRUE, col.names=FALSE)
