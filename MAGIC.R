pip install --user magic-impute
pip install --user phate


library(Rmagic)
library(ggplot2)
library(viridis)
library(phateR)

impute_seurat_MAGIC(seurat_object)

data(magic_testdata)
MAGIC_data <- magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"))
ggplot(MAGIC_data) +
  geom_point(aes(x=VIM, y=CDH1, color=ZEB1))