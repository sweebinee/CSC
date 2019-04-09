start on 25Mar
BC_Epithelial_cell_data + BC_Hematopoietic_cell_data

library(Seurat)
library(ggplot2)
library(RColorBrewer)

setwd('/storage2/Project/CSC/10X/DGIST_data02')
non_blood =readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Epithelial_cell_data/EpithelialCells_SeuratSet.rds")
blood = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/BC_blood_seurat.rds")
#> non_blood@meta.data$cell_label
#NULL
#> blood@meta.data$cell_label
#"BloodCell"


