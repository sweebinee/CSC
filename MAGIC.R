pip install --user magic-impute
pip install --user phate


library(Rmagic)
library(ggplot2)
library(viridis)
library(phateR)
library(scFunctions)

blood = readRDS(file="/storage2/Project/CSC/10X/DGIST_data02/BC_Hematopoietic_cell_data/BC_blood_seurat.rds")
impute_seurat_MAGIC(blood)

data(magic_testdata)
MAGIC_data <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
ggplot(MAGIC_data) +
  geom_point(aes(x=VIM, y=CDH1, color=ZEB1))