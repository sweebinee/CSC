# cut off
#source("http://cf.10xgenomics.com/supp/cell-exp/rkit-install-2.0.0.R")
setwd('/storage2/Project/CSC/10X/05_QC/PE2324252629')

library(cellrangerRkit) 
cellranger_pipestance_path = "/storage2/Project/CSC/10X/04_aggr/PE2324252629_aggr"

sample_lee <- load_cellranger_matrix_h5(cellranger_pipestance_path, barcode_filtered = FALSE)  # loading raw Gene Barcode Matrix

# plot_barcode_counts(gbm)
count_lee = exprs(sample_lee) #63677 x 2211840 sparse Matrix of class "dgTMatrix"
colSum_lee = colSums(count_lee)
rm(sample_lee)

colSumReduced = colSum_lee[colSum_lee>50] # 390 gives us ~50000 barcodes

# rank plots by condition
PE23_1 = colSumReduced[grepl("(^(\\w*)-1$)", names(colSumReduced))]
PE23_2 = colSumReduced[grepl("(^(\\w*)-2)", names(colSumReduced))]
PE24 = colSumReduced[grepl("(^(\\w*)-3$)", names(colSumReduced))]
PE25 = colSumReduced[grepl("(^(\\w*)-4$)", names(colSumReduced))]
PE26 = colSumReduced[grepl("(^(\\w*)-5$)", names(colSumReduced))]
PE29 = colSumReduced[grepl("(^(\\w*)-6$)", names(colSumReduced))]

png("PE23_1_50.png")
plot(1:length(PE23_1), sort(PE23_1, decreasing=T), pch=20, log="y")
dev.off()

png("PE23_2_50.png")
plot(1:length(PE23_2), sort(PE23_2, decreasing=T), pch=20, log="y")
dev.off()

png("PE24_50.png")
plot(1:length(PE24), sort(PE24, decreasing=T), pch=20, log="y")
dev.off()

png("PE25_50.png")
plot(1:length(PE25), sort(PE25, decreasing=T), pch=20, log="y")
dev.off()

png("PE26_50.png")
plot(1:length(PE26), sort(PE26, decreasing=T), pch=20, log="y")
dev.off()

png("PE29_50.png")
plot(1:length(PE29), sort(PE29, decreasing=T), pch=20, log="y")
dev.off()

elbow_point <- function(colsum, min, max){
  sortedGBM = sort(colsum, decreasing=T)
  point1 = c(min, as.integer(sortedGBM[min]))
  point2 = c(max, as.integer(sortedGBM[max]))
  
  lineVec = point2 - point1
  lineVecNorm = lineVec / sqrt(sum(lineVec^2))
  
  distances <- c()
  for(i in min:max){
    P = c(i, as.integer(sortedGBM[i]))
    pointVec = P - point1
    distVec = pointVec - (pointVec %*% lineVecNorm)*lineVecNorm
    distances <- c(distances, sum(distVec^2))
  }
  y_max_dist = sortedGBM[min + which.max(distances)]
  return(y_max_dist)
}

elbow_point(PE23_1,9000,12000) #238
elbow_point(PE23_2,15000,18000) #688
elbow_point(PE24,10000,14000) #509
elbow_point(PE25,9000,12000) #108
elbow_point(PE26,9000,12000) #106 
elbow_point(PE29,12000,20000) #239

rm(PE23_1)
rm(PE23_2)
rm(PE24)
rm(PE25)
rm(PE26)
rm(PE29)

# reconstruct countBasal and countCL  count matrix
colSumPE23_1 = colSum_lee[colSum_lee>238]
countFilteredPE23_1 = count_lee[,names(colSumPE23_1[colSumPE23_1>238])]
countFilteredPE23_1= countFilteredPE23_1[,grepl("(^(\\w*)-1$)", names(colSumPE23_1))]
countPE23_1 = as.matrix(countFilteredPE23_1)
rm(colSumPE23_1)
rm(countFilteredPE23_1)

colSumPE23_2 = colSum_lee[colSum_lee>688]
countFilteredPE23_2 = count_lee[,names(colSumPE23_2[colSumPE23_2>688])]
countFilteredPE23_2= countFilteredPE23_2[,grepl("(^(\\w*)-2$)", names(colSumPE23_2))]
countPE23_2 = as.matrix(countFilteredPE23_2)
rm(colSumPE23_2)
rm(countFilteredPE23_2)

colSumPE24 = colSum_lee[colSum_lee>509]
countFilteredPE24 = count_lee[,names(colSumPE24[colSumPE24>509])]
countFilteredPE24= countFilteredPE24[,grepl("(^(\\w*)-3$)", names(colSumPE24))]
countPE24 = as.matrix(countFilteredPE24)
rm(colSumPE24)
rm(countFilteredPE24)

countPE2324=cbind(countPE23_1, countPE23_2, countPE24)
rm(countPE23_1)
rm(countPE23_2)
rm(countPE24)

colSumPE25 = colSum_lee[colSum_lee>108]
countFilteredPE25 = count_lee[,names(colSumPE25[colSumPE25>108])]
countFilteredPE25= countFilteredPE25[,grepl("(^(\\w*)-4$)", names(colSumPE25))]
countPE25 = as.matrix(countFilteredPE25)
rm(colSumPE25)
rm(countFilteredPE25)

colSumPE26 = colSum_lee[colSum_lee>106]
countFilteredPE26 = count_lee[,names(colSumPE26[colSumPE26>106])]
countFilteredPE26= countFilteredPE26[,grepl("(^(\\w*)-5$)", names(colSumPE26))]
countPE26 = as.matrix(countFilteredPE26)
rm(colSumPE26)
rm(countFilteredPE26)

countPE2526=cbind(countPE25, countPE26)
rm(countPE25)
rm(countPE26)
countPE23242526=cbind(countPE2324, countPE2526)
rm(countPE2324)
rm(countPE2526)

colSumPE29 = colSum_lee[colSum_lee>239]
countFilteredPE29 = count_lee[,names(colSumPE29[colSumPE29>239])]
countFilteredPE29= countFilteredPE29[,grepl("(^(\\w*)-6$)", names(colSumPE29))]
countPE29 = as.matrix(countFilteredPE29)
rm(colSumPE29)
rm(countFilteredPE29)

rm(count_lee)
rm(colSum_lee)
rm(colSumReduced)

count_CSC = cbind(countPE23242526,countPE29)
rm(countPE29)

countFiltered = count_CSC
rm(count_CSC)
colnames(countFiltered) = gsub("^(\\w*)-1$", "PE23-1-\\1", colnames(countFiltered))
colnames(countFiltered) = gsub("^(\\w*)-2$", "PE23-2-\\1", colnames(countFiltered))
colnames(countFiltered) = gsub("^(\\w*)-3$", "PE24-\\1", colnames(countFiltered))
colnames(countFiltered) = gsub("^(\\w*)-4$", "PE25-\\1", colnames(countFiltered))
colnames(countFiltered) = gsub("^(\\w*)-5$", "PE26-\\1", colnames(countFiltered))
colnames(countFiltered) = gsub("^(\\w*)-6$", "PE29-\\1", colnames(countFiltered))

dim(countFiltered)
[1] 63677 73163

save(countFiltered, file = "CSC_countFiltered_method.RData")
write.csv(countFiltered, file = "count_CSC_subin_method.csv")

#####################################################
#####QC
##################################################
# C++ 11 compiler 필요 
# R >= 3.5 
source("http://bioconductor.org/workflows.R")
workflowInstall("simpleSingleCell")
# 안되면 biocLite() 해보자 + update all
source("https://bioconductor.org/biocLite.R")
biocLite("simpleSingleCell")
biocLite("scater")
biocLite("scran")
#`maximal number of DLLs reached... << 이런 에러 나면
# in R homedir ~ R_MAX_NUM_DLLS=256 (in MacOS)

library(scater)
library(biomaRt)
data_path = "/storage2/Project/CSC/10X/05_QC/PE2324252629"
figure_path = "/storage2/Project/CSC/10X/05_QC/PE2324252629"
setwd(data_path)
load(file = "CSC_countFiltered_method.RData")

countFiltered = countFiltered[rowSums(countFiltered) > 0, ]
dim(countFiltered)
[1] 36333 73163
cd = as.matrix(countFiltered)
logcd = log2(cd + 1)
dim(logcd)
[1] 36333 73163

cellInfo = data.frame(cell=colnames(countFiltered), 
                      condition=c(sapply(colnames(countFiltered)[grepl("PE23",colnames(countFiltered))], function(x) substr(x, 1, 6)),
                                  sapply(colnames(countFiltered)[grepl("PE24",colnames(countFiltered))], function(x) substr(x, 1, 4)),
                                  sapply(colnames(countFiltered)[grepl("PE25",colnames(countFiltered))], function(x) substr(x, 1, 4)),
                                  sapply(colnames(countFiltered)[grepl("PE26",colnames(countFiltered))], function(x) substr(x, 1, 4)),
                                  sapply(colnames(countFiltered)[grepl("PE29",colnames(countFiltered))], function(x) substr(x, 1, 4))
                                  )
)

BC_sce <- SingleCellExperiment(assays=list(counts=cd, logcounts=logcd), colData = cellInfo)
colnames(colData(BC_sce))
#[1] "cell"      "condition"

rm(countFiltered)
rm(cd)
rm(logcd)

#setwd(data_path)
#date = Sys.Date()
#ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="www.asia.ensembl.org")
#ensemblGenes <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name',  'chromosome_name', 'gene_biotype'), mart=ensembl)
#rownames(ensemblGenes) <- ensemblGenes[,1]
#mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]
#save(ensemblGenes, file=paste("ensemblGenes",date,".RData",sep=""))
load("ensemblGenes2018-01-14.RData")

mtGenes = ensemblGenes[ensemblGenes[,3]=="MT",]
is.mito = rownames(BC_sce) %in% mtGenes[,1]
length(is.mito[is.mito == TRUE])
#[1] 34
#rm(ensembl)
#rm(ensemblGenes)

BC_sce = calculateQCMetrics(BC_sce,  feature_controls=list(Mt=is.mito))
BC_sce = runPCA(BC_sce, pca_data_input = "colData")  
rm(cellInfo)

##[Before QC]#################################################################################
#[Before QC Plot]##
cellColor = c("#d73027","#d72776","#ffa500","#daff00","#1d6d3d","#4575b4","#800080") ##red,pink,orange,yellow,green,blue,purple
png("Before_QC_totalGeneCount_MTProp.png")
par(cex.axis=1.5, cex.lab=1.5)
plot(NULL, xaxt="n", xlab="Number of UMI counts",
     ylab="Percentage of reads mapped to mitochondrial genes", xlim=c(10^2,1*10^7), ylim=c(0,100), log="x")
axis( 1, 10^(seq(2, 7, 2)), c(expression(10^2), expression(10^4), expression(10^6)) )
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE23-1", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE23-1", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[1])
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE23-2", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE23-2", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[2])
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE24", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE24", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[3])
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE25", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE25", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[4])
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE26", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE26", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[5])
points(colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE29", "total_counts"],
       colData(BC_sce)[colData(BC_sce)[,"condition"]=="PE29", "pct_counts_Mt"], pch=20, cex=0.3, col=cellColor[6])
abline(h=15, v=c(1000,5000), col="red", lty=2, lwd=2) #h=10
legend(100000, 100,c("PE23-1","PE23-2","PE24","PE25","PE26","PE29"),
       pch=20,
       col=cellColor, bty="n", cex=1.5)
dev.off()

#[Before PCA Plot]##
png("Before_QC_PCA.png")
cellColor = c("#d73027","#d72776","#ffa500","#daff00","#1d6d3d","#4575b4","#800080")
PC1 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(BC_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(BC_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE23-1", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE23-1", "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE23-2", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE23-2", "PC2"], pch=20, col=cellColor[2], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE24", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE24", "PC2"], pch=20, col=cellColor[3], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE25", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE25", "PC2"], pch=20, col=cellColor[4], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE26", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE26", "PC2"], pch=20, col=cellColor[5], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE29", "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"condition"]=="PE29", "PC2"], pch=20, col=cellColor[6], cex=0.5)
legend(10, 14, c("PE23-1","PE23-2","PE24","PE25","PE26","PE29"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

#[MT-QC threshold decision Plot]##
pMT = 15

png("Before_QC_mtProp15_PCA.png")
cellColor = c("#d73027", "#4575b4") #red, blue
PC1 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(BC_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(BC_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_counts_Mt"]<pMT, "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_counts_Mt"]<pMT, "PC2"], pch=20, col=cellColor[2], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_counts_Mt"]>=pMT, "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_counts_Mt"]>=pMT, "PC2"], pch=20, col=cellColor[1], cex=0.5)
legend(0, 15, c(">15% of MT", "<15% of MT"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

# calculate pct_dropout
assay(BC_sce, "is_exprs") <- calcIsExprs(BC_sce,lowerDetectionLimit = 5)
colData(BC_sce)$pct_dropout <- (1 - apply(assay(BC_sce, "is_exprs"), 2, mean))* 100
colnames(colData(BC_sce))
[1] "cell"                                  
 [2] "condition"                             
 [3] "total_features"                        
 [4] "log10_total_features"                  
 [5] "total_counts"                          
 [6] "log10_total_counts"                    
 [7] "pct_counts_top_50_features"            
 [8] "pct_counts_top_100_features"           
 [9] "pct_counts_top_200_features"           
[10] "pct_counts_top_500_features"           
[11] "total_features_endogenous"             
[12] "log10_total_features_endogenous"       
[13] "total_counts_endogenous"               
[14] "log10_total_counts_endogenous"         
[15] "pct_counts_endogenous"                 
[16] "pct_counts_top_50_features_endogenous" 
[17] "pct_counts_top_100_features_endogenous"
[18] "pct_counts_top_200_features_endogenous"
[19] "pct_counts_top_500_features_endogenous"
[20] "total_features_feature_control"        
[21] "log10_total_features_feature_control"  
[22] "total_counts_feature_control"          
[23] "log10_total_counts_feature_control"    
[24] "pct_counts_feature_control"            
[25] "total_features_Mt"                     
[26] "log10_total_features_Mt"               
[27] "total_counts_Mt"                       
[28] "log10_total_counts_Mt"                 
[29] "pct_counts_Mt"                         
[30] "is_cell_control"                       
[31] "pct_dropout"  

#[Dropdout threshold decision Plot]##
pDROPOUT = 99.9
# table(rowData(BC_sce)$pct_dropout_counts < 99.0)

png("Before_QC_dropout_PCA_999.png")
cellColor = c("#d73027", "#4575b4")
PC1 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(BC_sce@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(BC_sce@reducedDims)$listData$PCA[,"PC1"], 
     attributes(BC_sce@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_dropout"]>=pDROPOUT, "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_dropout"]>=pDROPOUT, "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_dropout"]<pDROPOUT, "PC1"],
       attributes(BC_sce@reducedDims)$listData$PCA[colData(BC_sce)[,"pct_dropout"]<pDROPOUT, "PC2"], pch=20, col=cellColor[2], cex=0.5)
legend(-85, 15, c(">99.9% of dropout", "<99.9% of dropout"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

##[After QC]#################################################################################
#scesetFiltered = BC_sce[rowData(BC_sce)$pct_dropout_counts<99.9, colData(BC_sce)[,"pct_counts_Mt"]<pMT & colData(BC_sce)[,"pct_dropout"]<pDROPOUT]

##take a long time
#scesetFiltered_totcount1000 = BC_sce[rowData(BC_sce)$pct_dropout_counts<99.9, colData(BC_sce)[,"pct_counts_Mt"]<15 & colData(BC_sce)[,"total_counts"]> 1000]

#is.mito = rownames(scesetFiltered_totcount1000) %in% mtGenes[,1]
#scesetFiltered_totcount1000 = calculateQCMetrics(scesetFiltered_totcount1000,  feature_controls=list(Mt=is.mito))
#scesetFiltered_totcount1000 = runPCA(scesetFiltered_totcount1000, pca_data_input = "colData")  

scesetFiltered_totcount5000 = BC_sce[rowData(BC_sce)$pct_dropout_counts<99.9, colData(BC_sce)[,"pct_counts_Mt"]<15 & colData(BC_sce)[,"total_counts"]> 5000]

is.mito = rownames(scesetFiltered_totcount5000) %in% mtGenes[,1]
scesetFiltered_totcount5000 = calculateQCMetrics(scesetFiltered_totcount5000,  feature_controls=list(Mt=is.mito))
scesetFiltered_totcount5000 = runPCA(scesetFiltered_totcount5000, pca_data_input = "colData")  


dim(scesetFiltered_totcount5000)
[1] 19107 21557

scesetFiltered_PE23 = scesetFiltered[,colnames(scesetFiltered)[grepl("PE23",colnames(scesetFiltered))]]
saveRDS(scesetFiltered_PE23, "scesetFiltered5000_PE23.rds")
scesetFiltered_PE24 = scesetFiltered[,colnames(scesetFiltered)[grepl("PE24",colnames(scesetFiltered))]]
saveRDS(scesetFiltered_PE24, "scesetFiltered5000_PE24.rds")
scesetFiltered_PE25 = scesetFiltered[,colnames(scesetFiltered)[grepl("PE25",colnames(scesetFiltered))]]
saveRDS(scesetFiltered_PE25, "scesetFiltered5000_PE25.rds")
scesetFiltered_PE26 = scesetFiltered[,colnames(scesetFiltered)[grepl("PE26",colnames(scesetFiltered))]]
saveRDS(scesetFiltered_PE26, "scesetFiltered5000_PE26.rds")
scesetFiltered_PE29 = scesetFiltered[,colnames(scesetFiltered)[grepl("PE29",colnames(scesetFiltered))]]
saveRDS(scesetFiltered_PE29, "scesetFiltered5000_PE29.rds")

#saveRDS(scesetFiltered_totcount1000, "scsetFiltered1000_PCTDROPOUT999.rds")
saveRDS(scesetFiltered_totcount5000, "scsetFiltered5000_PCTDROPOUT999.rds")

#colData(BC_sce)

#logcd<-read.csv(file="logCSC.csv")
#scesetFiltered_totcount5000=readRDS(file="scsetFiltered5000_PCTDROPOUT999.rds")
#cell = scesetFiltered_totcount5000$cell 
#inferCNV <- subset(logcd, select=cell,sort=F)
#> ncol(inferCNV)
#[1] 11352
#> ncol(logcd)
#[1] 37384
#write.csv(inferCNV,file="inferCNV.csv")

#[After QC Plot]##
cellColor = c("#d73027","#d72776","#ffa500","#daff00","#1d6d3d","#4575b4","#800080") ##red,pink,orange,yellow,green,blue,purple
png("After_QC_totalGeneCount_MTProp.png")
par(cex.axis=1.5, cex.lab=1.5)
plot(NULL, xaxt="n", xlab="Number of UMI counts",
     ylab="Percentage of reads mapped to mitochondrial genes", xlim=c(10^2,1*10^7), ylim=c(0,100), log="x")
axis( 1, 10^(seq(2, 7, 2)), c(expression(10^2), expression(10^4), expression(10^6)) )
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-1", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-1", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[1])
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-2", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-2", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[2])
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE24", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE24", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[3])
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE25", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE25", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[4])
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE26", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE26", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[5])
points(colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE29", "total_counts"],
       colData(scesetFiltered_totcount5000)[colData(scesetFiltered_totcount5000)[,"condition"]=="PE29", "pct_counts_Mt"], pch=20, cex=1, col=cellColor[6])
abline(h=15, v=5000, col="red", lty=2, lwd=2)
legend(5000, 100,c("PE23-1","PE23-2","PE24","PE25","PE26","PE29"),
       pch=20,
       col=cellColor, bty="n", cex=1.5)
dev.off()

#[After PCA Plot]##
png("After_QC_PCA.png")
cellColor = c("#d73027","#d72776","#ffa500","#daff00","#1d6d3d","#4575b4","#800080")
PC1 <- format(attributes(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA)$percentVar[1]*100,digits=3)
PC2 <- format(attributes(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA)$percentVar[2]*100,digits=3)
par(mar=c(5,5,1,1), cex.axis=1.5, cex.lab=1.5)
plot(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[,"PC1"], 
     attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[,"PC2"], pch=20, col="black",
     xlab=paste("PC1 ", "(", PC1, "%)", sep=""), ylab=paste("PC2 ", "(", PC2, "%)", sep=""), cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-1", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-1", "PC2"], pch=20, col=cellColor[1], cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-2", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE23-2", "PC2"], pch=20, col=cellColor[2], cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE24", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE24", "PC2"], pch=20, col=cellColor[3], cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE25", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE25", "PC2"], pch=20, col=cellColor[4], cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE26", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE26", "PC2"], pch=20, col=cellColor[5], cex=0.5)
points(attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE29", "PC1"],
       attributes(scesetFiltered_totcount5000@reducedDims)$listData$PCA[colData(scesetFiltered_totcount5000)[,"condition"]=="PE29", "PC2"], pch=20, col=cellColor[6], cex=0.5)
legend(50, 20 , c("PE23-1","PE23-2","PE24","PE25","PE26","PE29"),
       pch=20,
       col=cellColor, bty="n", cex=1)
dev.off()

############
library(scater)
library(scran)

# scesetFiltered=readRDS(file="scsetFiltered1000_PCTDROPOUT999.rds")
#scesetFiltered=readRDS(file="scesetFiltered5000_PE23.rds")
#scesetFiltered=readRDS(file="scesetFiltered5000_PE24.rds")
#scesetFiltered=readRDS(file="scesetFiltered5000_PE25.rds")
#scesetFiltered=readRDS(file="scesetFiltered5000_PE26.rds")
#scesetFiltered=readRDS(file="scesetFiltered5000_PE29.rds")

#take very very long time
clusters <- quickCluster(scesetFiltered)
# 19107 genes 
# save(clusters, file ="scran_clusters_1000.RData")
#save(clusters, file ="scran_clusters_5000_PE23.RData")
#save(clusters, file ="scran_clusters_5000_PE24.RData") #4385 cells / 5clusters
#save(clusters, file ="scran_clusters_5000_PE25.RData") #2015 cells / 3clusters
#save(clusters, file ="scran_clusters_5000_PE26.RData") #2705 cells / 5clusters
#save(clusters, file ="scran_clusters_5000_PE29.RData") #5485 cells / 4clusters

#clusters <- load("scran_clusters_5000_PE26.RData")

# scesetFiltered <- computeSumFactors(scesetFiltered, cluster=clusters, sizes=seq(500,600,10))
scesetFiltered <- computeSumFactors(scesetFiltered, cluster=clusters)
#computeSumFactors()
#Methods to normalize single-cell RNA-seq data by deconvolving size factors from cell pools.
scesetFiltered <- normalize(scesetFiltered)

#pdf("Figure_SF_LibrarySize_1000.pdf")
png("PE26_SF_LibrarySize_5000.png")
plot(main="PE26",sizeFactors(scesetFiltered), scesetFiltered$total_counts/5000, log="xy",
     ylab="Library size (kilos)", xlab="Size factor", pch=20)
dev.off()

#saveRDS(scesetFiltered, "scsetFiltered1000_PCTDROPOUT999_norm.rds")
saveRDS(scesetFiltered, "PE26_scsetFiltered5000_PCTDROPOUT999_norm.rds")
#scesetFiltered = readRDS("PE26_scsetFiltered5000_PCTDROPOUT999_norm.rds")

library("aroma.light")
var.fit <- trendVar(scesetFiltered, method="spline", parametric=TRUE, use.spikes=FALSE, span=0.2)#, design = batchDesign)
var.out <- decomposeVar(scesetFiltered, var.fit)
hvg <- var.out[which(var.out$FDR <= 0.05 & var.out$bio > .5),]
dim(hvg)
# PE24: 696   
# PE25: 368   
# PE26: 372
# PE29: 591 

#pdf("Figure_MeanExp_VarExp_BIO_0.5_1000.pdf")
png("PE29_MeanExp_VarExp_BIO_0.1_5000.png")
plot(main="PE29",var.out$mean, var.out$total, pch=16, cex=0.3, xlab="Mean log-expression",ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
points(var.out$mean[var.out$FDR <= 0.05 & var.out$bio > .1], var.out$total[var.out$FDR <= 0.05 & var.out$bio > .1], pch=16, cex=0.3, col="red")
dev.off()

save(hvg, file="PE26_hvg_5000.RData")


#####################################################
#####codeSeurat
#####################################################
source("http://bioconductor.org/biocLite.R")
biocLite("Seurat")

library(Seurat)
library(scater)
library(dplyr)

#load("hvg_1000.RData")
load("PE24_hvg_5000.RData")

#scesetFiltered = readRDS(file = "scsetFiltered1000_PCTDROPOUT999_norm.rds")
scesetFiltered = readRDS(file = "PE24_scsetFiltered5000_PCTDROPOUT999_norm.rds")

sce2seurat <- setClass("sce2seurat",
  slots=c(
    data= "ANY",
    meta.data = "data.frame",
    scale.data="ANY",
    cell.names = "vector",
    var.genes="vector",
    is.expr = "numeric",
    ident="factor",
    snn= "dgCMatrix",
    dr = "list",
    calc.params="list"
  )
)

#original code wasn't working
#BC_seurat = sce2seurat(data = Matrix(exprs(scesetFiltered), sparse = T),
#                            scale.data = exprs(scesetFiltered),
#                            cell.names = colnames(exprs(scesetFiltered)), 
#                            var.genes = rownames(scesetFiltered[rownames(hvg),]))

#exprs(object): returns the matrix of (log-counts) expression values
#https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/scater/inst/doc/vignette-transition.html
BC_seurat = sce2seurat(data = Matrix(logcounts(scesetFiltered), sparse = T),
                            scale.data = logcounts(scesetFiltered),
                            cell.names = colnames(scesetFiltered), 
                            var.genes = rownames(scesetFiltered[rownames(hvg),]))
num.genes <- colSums(BC_seurat@scale.data > 0)
num.mol <- colSums(BC_seurat@scale.data)
cells.use <- names(nuSeurat.genes[which(num.genes > 0)])
nGene <- num.genes[cells.use]
nUMI <- num.mol[cells.use]
BC_seurat@meta.data <- data.frame(nGene, nUMI)

##### Dimension reduction

PCA = 50
BC_seurat <- RunPCA(BC_seurat, pcs.compute = PCA, weight.by.var = FALSE)

BC_seurat <- RunTSNE(BC_seurat, dims.use = 1:PCA, do.fast = T, seed.use = 42, perplexity=20)
BC_seurat <- FindClusters(BC_seurat, reduction.type="pca", dims.use = 1:PCA, save.SNN = TRUE, force.recalc = TRUE)

pdf("PE29_seurat_pca.pdf")
PCAPlot(BC_seurat)
dev.off()

pdf("PE29_seurat_tsne.pdf")
TSNEPlot(BC_seurat)
dev.off()

save(BC_seurat, file="PE24_seurat_cluster_5000.RData")
# save(BC_seurat, file="seurat_cluster_1000.RData")

#check sdev
# load(file="seurat_cluster_1000.RData")
# load(file="seurat_cluster_5000.RData")

#BC_seurat@dr$pca@cell.embeddings
pdf("PE29_PCA_Seurat_cluster_5000.pdf") 
df = data.frame(x=BC_seurat@dr$pca@cell.embeddings[, "PC1"], 
                y=BC_seurat@dr$pca@cell.embeddings[, "PC2"], 
                expression=BC_seurat@ident)
centers = data.frame(df %>% dplyr::group_by(BC_seurat@ident) %>% summarize(x = median(x = x), 
                                                                                y = median(x = y)))
ggplot() + 
  geom_point(data=df, size=1, aes(x=x, y=y, colour=expression)) + 
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) +
  geom_text(data = centers,
            mapping = aes(x = x, y = y, label = BC_seurat.ident),
            size = 5)
dev.off()


pdf("PE29_tSNE_Seurat_cluster_5000.pdf") 
df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=BC_seurat@ident)
centers = data.frame(df %>% dplyr::group_by(BC_seurat@ident) %>% summarize(x = median(x = x), 
                                                                                y = median(x = y)))
ggplot() + 
  geom_point(data=df, size=0.5, aes(x=x, y=y, colour=expression)) + 
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) +
  geom_text(data = centers,
            mapping = aes(x = x, y = y, label = BC_seurat.ident),
            size = 5)
dev.off()

################################################################################################
# inferCNV 로 찾아낸 cancer population plot그리기 생략
################################################################################################
setwd('/storage2/Project/CSC/10X/05_QC/PE2324252629')

load("PE25_seurat_cluster_5000.RData")

#it takes long time 10 min
BC_seurat.markers <- FindAllMarkers(object = BC_seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25, min.cells.gene=2)

#save(BC_seurat.markers, file="PE24_BC_5000_Markers.RData")
load(file="PE25_BC_5000_Markers.RData")
load("ensemblGenes2018-01-14.RData")

BC_seurat.markers.genesymbol = merge(BC_seurat.markers, ensemblGenes, by.x="gene", by.y="ensembl_gene_id", sort=FALSE, all.x=TRUE)
BC_seurat.markers.genesymbol = BC_seurat.markers.genesymbol[BC_seurat.markers.genesymbol["p_val_adj"] < 0.01,]
write.csv(BC_seurat.markers.genesymbol, file="PE26_BC_5000_seurat.markers.genesymbol.csv")

> colnames(BC_seurat.markers.genesymbol)
 [1] "gene"               "p_val"              "avg_logFC"         
 [4] "pct.1"              "pct.2"              "p_val_adj"         
 [7] "cluster"            "external_gene_name" "chromosome_name"   
[10] "gene_biotype"

######################################################################################
scesetFiltered = readRDS(file = "PE25_scsetFiltered5000_PCTDROPOUT999_norm.rds")

colData(scesetFiltered)$cluster = BC_seurat@ident

table(colData(scesetFiltered)[,"cluster"]==4)
length(which(colData(scesetFiltered)[,"cluster"]==0))

library(RColorBrewer)

################################### gene plot
# ERBB2 EGFR BRCA1 BRCA2 PTEN TP53 ESR1 RB1 PTPRC
#CD8+ naive T-cell
TGgene = c('AMBN','BCAS2','BTNL8','BUD31','C19orf53','CA14','CCDC87','CCR8','CD8A','CD8B','CDK5RAP1','CGRRF1','CILP','COX4I1','CWF19L1','DDX24','EEF1D','FXYD7','GDAP2','GJB4','GPR15','GPR52','HAUS3','HIGD2A','HIST1H3A','HIST1H4F','HTR1B','IL21R','JMJD6','JOSD1','KIAA1109','KRT1','LIN28A','LUZP4','MOGAT2','MS4A5','MYL1','NDUFA1','NDUFA4','NDUFS5','NGDN','NKTR','OMG','PCIF1','PSG11','RFX2','RNF7','RNMT','RPP38','RRH','RRP8','SETD2','SHFM1','SKI','SLC17A4','SLC35E1','SMCP','SMR3B','SON','SS18L2','TNKS2','TPP2','USP36','ZNHIT3')
#CD8+ T-cell
TGgene = c("AAK1","APBB1","ARHGEF1","BTN2A1","C7orf26","CA6","CASP8","CBY1","CCDC25","CCDC53","CCR7","CD160","CD27","CD3D","CD7","CD8A","CD8B","CD96","CEPT1","CIAPIN1","CLUAP1","COG2","COPZ1","CRTAM","CTSW","CX3CR1","DHX15","DIDO1","DNAJB1","DPP8","DSC1","FAM134C","FBXW4","FKTN","FNBP4","FTO","GGNBP2","GIMAP4","GJC2","KLRB1","KLRG1","KRT2","LAIR2","LSM14A","LY9","MED17","MKRN2","MMP19","MSL3","MTRF1","MYOM1","NAA16","NDUFS2","NDFIP1","NFKB1","RBM34","RNF113A","RPL37A","RING1","RWDD3","S100B","SIRPG","SLC1A7","SSTR3","RBL2","SDAD1","SHANK1","TBCC","SFPQ","SDCCAG3","TMEM41B","TTN","TOMM7","TRAF3IP3","TSPAN32","ZBTB11","ZC3HAV1","ZNF154","ZNF611","ZNF639","ZNF200","UBE2Q1","UBQLN2","UTP20","WDR82","USP47","YLPM1","PRL","PRMT2","NPRL2","PCNT","PTPN4","PURA","PTGDR","PSD","RAPGEF6","RASA2","PRPF4B","PLCG1","PLXDC1","PFN2","POLR3E","NPAT","POP5","NKRF","IPCEF1","IRF3","GZMH","GZMK","HNRNPL","KLHL3","IL16","GZMM","HNRNPA0")
#CD8+ Tcm
TGgene = c('ABCD2','ADAT1','ADCYAP1R1','AKAP3','ALK','ATF7IP','ATXN7','BMPR1A','C21orf2','C21orf59','CASP8','CD27','CD28','CD3E','CD48','CD7','CD8A','CD8B','CD96','CDC14A','CGRRF1','CHST5','CRTAM','CTSW','CXCR3','CYP20A1','DCTN6','DNAJB1','DUSP11','DYNLT1','ELP3','ESR2','GIMAP4','GIMAP6','GNLY','GPR171','GZMA','GZMH','GZMK','HAO2','HAUS3','HELZ','HIST1H4F','HMGN4','IL2RB','INPP4A','INTS5','ISCA1','KIF2A','KLRD1','KLRG1','KRT1','LEPROTL1','LRIG2','MED6','MMP11','NAA16','NCK1','LAG3','PARP11','PCDHA10','PDCD1','PTPN4','PPWD1','PTPRCAP','RASA2','RPS6KB1','S100B','RGS1','SEC24A','SFXN1','SGCD','SLC6A7','SMAP1','STUB1','SIT1','SH2D1A','SYN2','TBCC','TMEM30B','TNFSF8','TNKS2','TNRC6B','TPP2','TRADD','TRAF3IP1','VPS33A','WDR18','XCL1','USP36','ZAP70','ZBTB1','ZC3H13','ZMYND11')
#CD8+ Tem
TGgene = c('AAK1','ABCD2','ABCF1','ABCF3','ABT1','ACADVL','ACTR1B','ADAT1','AHNAK','ALKBH4','AMBRA1','ANAPC2','ANGEL2','ARFRP1','ARHGEF1','ARHGEF5','ARPC5L','ASTE1','ATR','B4GALT3','BIN2','BIN3','BMPR1A','BTN2A1','C14orf169','C1orf35','C21orf2','C2orf42','C7orf26','C8G','CACNB1','CALM1','CAPN10','CAPN2','CAPZB','CASP8','CCDC130','CCDC85C','CCL4','CCR5','CD160','CD2','CD300A','CD3D','CD3E','CD3G','CD7','CD8A','CD8B','CD96','CDC14A','CDK10','CDKN2AIP','CEP250','CHD8','CHST12','CIAO1','CNNM2','COL5A3','COLQ','COPB1','CORO7','CREBZF','CROCC','CRTAM','CSNK1G2','CSTF2T','CTBP1','CTDSP1','CTR9','CTSW','CX3CR1','CXCR3','CXCR6','CYTH4','DAXX','DDX18','DEFB126','DHX16','DHX8','DIAPH1','DIDO1','DLEC1','DMWD','DNAJB1','DNAJC24','DPP8','DUSP8','DYNC1H1','DYNLT1','E4F1','ELK4','ELP3','EMD','ERAL1','ERN1','EXOC2','FASLG','FBXO3','FBXO31','FBXW4','FLT4','FYCO1','GCC1','GIMAP4','GIMAP6','GIPR','GNLY','GOLGA1','GOLGA4','GPKOW','GPR171','GPR52','GPR65','GPR68','GRAP2','GSDMD','GTF3C1','GYG1','GZMA','GZMB','GZMH','GZMK','GZMM','HAUS3','HIVEP3','HLA-A','HLCS','HMGN4','HMGXB3','HMOX2','IDH3B','IFNG','IKBKAP','IKZF3','IL10RA','IL12RB1','IL18RAP','IL2RB','IMP3','INPP5E','IRF3','ITGAL','JMJD6','KIAA0196','KIAA1109','KIF22','KLHDC4','KLHL11','KLRB1','KLRD1','KLRG1','KRI1','LAG3','LAIR2','LIPT1','LTB4R2','LTK','LZTR1','MAML1','MAN2C1','MAP3K10','MAP7D1','MAPK13','MARK3','MDN1','MKNK1','MRFAP1L1','MRPL22','MRPS22','MSH3','MT2A','MTO1','MUS81','MYO1F','N4BP1','NCR1','NCR3','NDUFB1','NDUFB2','NDUFS6','NKG7','NMUR1','NOL6','NPRL2','NRBP1','OSBPL7','OTOF','OTUD7B','PABPC3','PANK4','PAPOLG','PCNT','PDCD1','PHKG2','PLCG1','PLXDC1','PMS1','PNMA3','POLG','PPP1CA','PPP2R5C','PRDM8','PREB','PRF1','PRKG2','PRMT7','PSMC5','PSME1','PSTPIP1','PTGDR','PTPN4','PTPRA','PTPRCAP','PUF60','PURA','PVRIG','PWWP2A','PYHIN1','PZP','RBL2','RFX1','RGS1','RGS9','RHOG','RHOT2','RIC8A','RIPK1','RNF113A','RNF126','RNF167','RNF25','RPN2','RPUSD2','RWDD1','S100B','SACM1L','SART3','SASH3','SBF1','SCAP','SEC16A','SEC24C','SF3A1','SF3B2','SGSM2','SH2D1A','SIT1','SLAMF1','SLC1A7','SLC25A12','SLC25A20','SLC35B1','SLC38A1','SNRNP200','SNTB2','SPSB3','SRCAP','SRRT','SSR2','STK25','STUB1','STUB1','STX18','STX4','TAF10','TBCC','TBX21','TCOF1','TIMM22','TMEM184B','TNIP2','TRAF3IP3','TRAPPC4','TRMT2A','TRMU','TRNAU1AP','TTC22','TUFM','TULP4','UBE2Q1','UBN1','UBTF','URB2','URGCP','USP1','USP34','USP47','VPS4A','WAS','WDR82','WSB2','WWP2','XCL1','ZAP70','ZBTB1','ZBTB39','ZMYM1','ZNF142','ZNF276','ZNF335','ZNF394','ZNF428','ZNF549','ZNF639','ZNF696','ZNF79')
#cDC
TGgene = c('ACTR3','ALCAM','ALDH1A2','ALOX15','ANXA1','CCDC88A','CCL13','CCL17','CCL23','CCL24','CD163','CD1A','CD1B','CD1C','CD1E','CD209','CD80','CD86','CD93','CLEC10A','CLEC4A','CRH','DBI','DNASE1L3','FCER1A','FGL2','FLT3','GFRA2','ITGAX','KCNK13','PITPNA','RAB7A','RRP1B','S100A10','SLAMF8','SSR1','TCTN3','WDFY3')
#Chondrocytes
TGgene = c('10-Sep','ABCC9','ABCD1','ABI2','ACAN','ACBD3','ADM2','ADRA1D','AK1','ANGPTL7','APOC3','ARCN1','ARL1','ARNT','BAG3','BCL2L13','C11orf95','CABP1','CACNA1C','CACNB1','CALCOCO2','CAMLG','CCL27','CCNB1','CDIPT','CDK2AP1','CELA2B','CILP','CLEC3B','CMPK1','CNIH4','COL10A1','COL14A1','COL5A3','COMP','COPA','COPB1','COPS8','CORIN','CRTAC1','CSNK1A1','CTRB2','CUL7','CYP19A1','DMWD','DPT','DPY19L4','DSTN','DVL1','DYNLRB1','EID1','ELN','EMILIN1','ENO1','ERC1','ERF','ERG','ETF1','FBXW11','FGF7','FNDC3A','FOXC2','FOXD2','FSHB','FTO','FZD9','GANAB','GDF5','GLG1','GLT8D1','GMPPA','GOLGA4','GORASP2','GOSR2','GRIA3','HAS1','HDAC11','HDLBP','HOXD3','HSPB7','IBSP','IDUA','IFNB1','IFT46','IGF1','IL17B','IL17RC','IQSEC2','IRGC','ISLR','KDELR1','KIAA0368','KLF15','KLHL26','KRT10','LAPTM4A','LBP','LEP','LPAR4','LRRC15','LTBR','MAP3K6','MAP4K5','MAST2','MCHR1','MIER2','MIF','MIOS','MKRN2','MLLT1','MLN','MORF4L2','MPHOSPH6','MYH13','MYOC','NAALADL1','NCKAP1','NCOR2','NDFIP1','NEUROG2','NFATC4','NFKBIL1','NGRN','NKX3-1','NPHP4','NPLOC4','NRBF2','NXPH3','OCRL','OGN','OLAH','OMD','OS9','OTOF','OTUD7B','P4HB','PCDHGC3','PDE3A','PFN2','PGK1','PGM1','PHLDB1','PJA2','PKNOX2','PODNL1','POFUT2','POMT2','PPIL6','PRDX4','PRELP','PRG4','PRKRA','PRL','PRRC1','PTH1R','PTP4A2','PTPN11','PYY','RAB3A','RANBP9','REM1','RGR','RGS11','RIC8A','RNF11','ROS1','RUNX1','RUNX2','S100A6','SCFD1','SCRG1','SEC61A1','SERINC3','SGCD','SLC13A4','SLC26A10','SLC38A3','SLC41A3','SLC5A12','SMAD5','SNED1','SNTB2','SNUPN','SNX13','SPAG9','SPATA7','SPIN1','SSX1','STAT2','SVEP1','TACR3','TBC1D16','TBC1D9B','TBX4','TCEAL4','TF','TM2D1','TM9SF1','TMED10','TMED7','TMEM59','TMEM59L','TNFSF11','TNNT3','TNXB','TRIM3','TSSK1B','TTC23','TTC37','TUB','TUBG2','UBXN6','WASL','WDR41','WISP1','WWP2','YAP1','YIF1A','YIPF2','ZBTB16','ZFAND3','ZFPL1','ZMYND11','ZNF358','ZNF471')
#Class-switched memory B-cells
TGgene = c('ABI1','ADAMDEC1','AFTPH','BAIAP3','CD80','COL19A1','CPSF4','CR1','CXCR5','DEPDC5','EPS15','FCRL2','GPR25','HLA-DQB2','KHDRBS2','MS4A1','NARFL','NDUFA9','NGLY1','PAX5','PIKFYVE','PRDM10','PTPN6','RAD17','RAPGEF1','SCRN1','SEC24A','SLC12A3','SNED1','SP140','SPIB','SUB1','TAF6','TNFRSF13B','TNFRSF17','TRAF3','TRRAP','UBE2G1','UBE2I','UBE2N','ZBTB32')
#aDC
TGgene = c('ABCA6','ABI1','ABTB2','ACHE','ACOT9','ADPGK','ADPRH','ALDH1A2','ALOX15B','ANKFY1','ANXA5','ARF3','ARFRP1','ARHGEF11','ARL8B','ARPC4','ATP1B3','AZIN1','BCKDK','BCL2L11','BCL2L13','BLVRA','C1orf27','C1QA','C1QB','C3AR1','C5orf15','CAMK1G','CASP5','CCL1','CCL13','CCL17','CCL18','CCL19','CCL22','CCL23','CCL24','CCL4','CCL7','CCL8','CCR5','CD209','CD300C','CD80','CD86','CHFR','CHMP5','CLIP1','CMKLR1','CUL1','CXCL13','CXCL9','CYTH4','DENND1A','DNAJA1','DNPEP','DOT1L','DYNLT1','EIF5','ELL','ENO1','ETV3','EXOC5','FAM175B','FBXL4','FCER1G','FCER2','FPR3','GMFB','GNG5','GPN2','GPR107','GPR132','GRB2','H6PD','HAMP','HCK','HLA-DQA1','HPS5','HRH2','HS3ST3B1','IL10','IL10RA','IL12B','IL19','IL2RA','IL3RA','IL9','IMPDH1','IRF4','KCNK13','KCNMB1','KDM6B','LAIR1','LILRA5','LILRB1','LILRB2','LILRB5','LOR','MAGT1','MAP2K1','MAP3K13','MED25','MFN1','MGRN1','MIIP','MMP25','MTF1','MTHFD2','N4BP1','NECAP2','NETO2','NEU3','NFE2L2','NFKB1','NFKBIB','NRAS','NSUN3','NUP62','OGFR','OPA3','OSM','P2RX7','PAK2','PGK1','PITPNA','PTGIR','RAB21','RAB35','RAB5A','RAB8A','RAMP3','RAPGEF1','RELA','RHOG','RIN2','RPGRIP1','RRP1','S100A10','SCARF1','SCYL2','SEC24A','SIGLEC1','SIGLEC7','SIGLEC9','SLAMF1','SLC1A2','SLC25A28','SLC6A12','SLCO5A1','SNX11','SOCS3','SPAG9','SRC','STAT2','TAX1BP1','TBC1D13','TBC1D22B','TCF21','TDRD7','TFEC','TMCC2','TMEM41B','TMSB10','TMX1','TNFRSF4','TNIP2','TOR1B','TPI1','TRAF1','TREX1','TRIP4','TRPC4AP','TXN','UBE2Z','VRK2','WTAP','XIAP','XPNPEP1','ZBTB17','ZNF654')
#Adipocytes
TGgene = c('ADH1B','ADIPOQ','ATP1A2','ATP5G3','C6','CILP','COL5A3','DBI','DLAT','ECHDC1','GPD1','HADHA','HP','LBP','PLIN1','PNPLA2','PPP1R1A','PPP2R1B','PTGER3','SLC25A6','TF')
#Astrocytes
TGgene = c('ACSS3','ACTA2','ACTG2','ADAM12','ADAM19','ADAMTS2','AFF3','AGPAT4','ANKRD1','APLP1','ASAP3','ATF7IP2','BGN','BST1','BTN3A3','C1S','C9orf16','CACNA1H','CAND2','CBR3','CBX2','CCL2','CDH11','CDH6','CENPI','CFH','CFI','CHN1','CHRDL1','CLCN2','CLDN11','CLSTN2','CLU','CNN1','COL11A1','COL13A1','COL1A1','COL1A2','COL3A1','COL5A1','COL5A2','COL6A1','COL6A2','COL6A3','COL8A2','COPZ2','CPA4','CRABP2','CREB3L1','CRLF1','CXCL12','DACT1','DAPK2','DCHS1','DCN','DHRS3','DIO2','DIRAS3','DNA2','DOK5','DPF3','DPYSL3','DUSP2','EDIL3','EFEMP2','FABP7','FAM64A','FAP','FBLN2','FBN1','FBXL7','FGFR1','FGG','FILIP1L','FKBP10','FLNC','FN1','FSD1','FZD2','FZD7','GABRQ','GATA6','GBP2','GEM','GFPT2','GINS4','GLI2','GLIPR1','GLT8D2','GNG11','GNG3','GPC4','GPX7','GREM1','GRM7','H2AFV','HEPH','HIST1H1D','HOXB2','HS3ST3A1','IGFBP2','IGFBP3','IGFBP5','IGFBP6','ITGA7','KIF18A','KIF20B','KIFC1','LAMA4','LIF','LOX','LOXL1','LRP4','LRRC17','LUM','LYPD1','MAP1A','MATN2','MCM3','MDK','MEG3','MEST','MFAP4','MFAP5','MGP','MMP2','MN1','MOXD1','MRC2','MXD3','MXRA5','MXRA8','MYL9','MYLK','NAP1L3','NCAPG','NES','NFASC','NID2','NNAT','NR2F1','NT5DC2','NT5DC3','OLFML2A','OLFML3','OXTR','PADI2','PALM','PAX6','PCOLCE','PCOLCE2','PDE5A','PDGFRA','PDGFRB','PLN','PNMA2','PNMAL1','POSTN','POU3F3','PPP1R3C','PSRC1','PTGIS','PTN','PTPRN','PTX3','RARB','RARRES1','RARRES2','RBP1','RCN3','REC8','RECQL4','RGS4','RGS5','RIBC2','SCARA3','SCN5A','SDC2','SEMA3E','SERPINF1','SERPINH1','SFRP4','SLC12A8','SLC14A1','SLC22A17','SPARC','SPON2','SRPX2','ST8SIA2','STAC','SULF1','SYNPO','SYT11','SYTL2','TAGLN','TFAP2C','TFPI','TGM2','THY1','TIMP1','TIMP2','TMEM135','TNC','TNFSF4','TNS1','TPM2','TRIM45','TRIOBP','TRO','TRPV2','UCHL1','VCAM1','VCAN','WNT5A','XRCC2','XYLT1','ZCCHC24','ZNF239')
#Bcell
TGgene = c('ACTN2','AFF2','AFTPH','AICDA','ANKMY1','AP3B1','ARHGAP17','ATF7IP','BAIAP3','BCL2L11','BLK','BMP8B','BTK','C12orf49','C5orf15','CCR6','CD180','CD19','CD22','CD37','CD53','CD72','CD79A','CD79B','CDC40','CEACAM21','CEPT1','CHAD','CIITA','CNOT1','CNR1','CNR2','COL19A1','CR1','CSNK1G3','CXCR5','DAXX','DCLRE1C','DEF8','DEPDC5','DNASE1','EGOT','FCRL2','GDI2','GGA2','GPR18','GPR25','GPRC5D','HDAC7','HLA-DOA','HLA-DPB1','HSPA4','HTR3A','IFNA2','IFNW1','IKZF3','IL17A','INPP5B','ITSN2','JMJD1C','KCNIP2','KCNN3','KIAA0430','KIAA1033','LSM6','LY86','LY9','MAP3K9','MBD4','MCM9','MFN1','MGAT5','MIOS','MMP17','MRM1','MS4A1','MYOT','NSUN5','NUP160','P2RY10','PAX5','PGR','PHKB','PIKFYVE','PLA2G2D','PNOC','PNOC','POLR3K','POU2F1','POU2F2','PRDM2','PRDM4','PRKCB','PWP1','QRSL1','RECQL5','RIC3','RNGTT','RPS11','RPS16','RRAS2','S1PR2','S1PR4','SEC62','SIPA1L3','SLC13A2','SLC24A1','SLC30A4','SMC6','SNX2','SP140','SPIB','STAG3','STAP1','SYPL1','TCL1A','TCL1B','TCL6','TERT','TLR7','TNFRSF13B','TNFRSF17','TRAF3','TROVE2','UBE2G1','UTP6','VPREB3','WDR11','ZNF154','ZNF202','ZNF208','ZNF37A','ZNF638','ZNF688','ZNF701','ZZZ3')
#Memory Bcell
TGgene = c('ACRV1','ADAM20','ADAM21','ADAM30','ADAMTS12','ADCY2','AICDA','AIPL1','ANKRD34C','AP1M2','AQP8','ART1','BAIAP3','BLK','C4BPA','CAPN3','CASQ2','CCR6','CCR9','CD180','CD19','CD1C','CD22','CD37','CD3EAP','CD72','CD79A','CD79B','CER1','CETP','CHP2','CHRM2','CHRNA2','CHST5','CLDN17','CNR2','COLEC10','COX6A2','CPA2','CPB1','CSN1S1','CXCL13','CXCR5','CYLC2','CYP2A7','CYP2C19','DCC','DDX4','DPP6','DSCR4','DSP','FCRL2','FMO1','FMO6P','FSCN2','FSCN3','FSHR','GABRA4','GAD2','GGA2','GK2','GLYAT','GNAT2','GNRHR','GPRC5D','GPX5','GRIN2B','GRM6','HCRTR2','HECW1','HLA-DPB1','HNRNPL','HRH4','HSD3B2','HTN3','HTR3A','INHBC','KALRN','KCNA5','KCNJ10','KHDRBS2','KIAA0125','KIF5A','KRT2','KRT75','LECT2','LY86','MAP3K9','MBD4','MBL2','MC4R','MEFV','MEP1B','MGAT5','MIOS','MOGAT2','MS4A1','MS4A5','MYOZ3','NCAN','NPHS2','NPY5R','NT5C','NTRK3','OBSCN','ODC1','OTC','PAX5','PGLYRP4','PIKFYVE','PLIN1','PNLIPRP1','PNOC','POU4F2','PRKCB','PROZ','PRPH2','PTH1R','QRSL1','RIC3','RNGTT','RRH','S100G','SCRN1','SELP','SERPINA4','SHISA6','SIGLEC6','SLC12A3','SLC17A1','SLC17A7','SLC24A2','SLC30A10','SLC5A7','SLCO1C1','SLN','SP140','SPIB','SSX3','STAP1','SYN2','SYPL1','TAS2R14','TCTN2','TMPRSS11D','TNFRSF13B','TNFRSF17','TRMT61A','TRPM3','TSHB','TSPAN13','ULK4','UNC5C','VN1R1','VPREB3','WNT16','WNT2','ZBTB32','ZIC3','ZNF548','ZNF747')
#Native Bcell
TGgene = c('08-Mar','ADAM20','AKAP6','AP3B1','BCL2L10','BCL2L11','BLK','BMP3','C10orf76','CACNA1F','CAPN3','CCR6','CD180','CD19','CD1A','CD22','CD37','CD72','CD79A','CD79B','CDK13','CIITA','COL19A1','CSNK1G3','CUBN','CXCR5','DAZL','DEF8','DSP','EGOT','FCER2','FCRL2','FRS2','GCM1','GGA2','GH1','GMFB','GNG3','GPR18','HLA-DOA','HSPA4','KHDRBS2','LY9','MAP3K9','MATN1','MBD4','MCM9','MFN1','MGAT5','MMP17','MS4A1','MYBPC2','MYO3A','N4BP3','NOC3L','P2RY10','PAX5','PGAM2','PHKG1','PIKFYVE','PNOC','POU2F1','PRDM2','PRDM4','PRKCB','PTCH2','PWP1','PYGM','RB1','RBM15','RERE','RRAS2','SDK2','SIPA1L3','SLC30A4','SMC6','SNTG2','SNX2','SP140','SPIB','STAG3','STAP1','SYN3','TBC1D5','TCL1A','TCL1B','TCL6','TRA2B','TRAPPC9','TREML2','TSPAN13','UBE2O','USP6','USP7','UTP6','VPREB3','WDR74','ZNF154')
#Pro-Bcell
TGgene = c('ADAM21','ADAMTS8','ADARB2','ADCY8','AFM','AGXT','AHDC1','AHSP','AKAP8L','ALOX15B','ALPL','ANAPC2','ANXA3','APOC3','ARG1','ARPP21','ART4','ASPM','ATP1B4','AZU1','BLK','BMP10','BMX','BTBD7','C16orf59','C2orf49','CA1','CA14','CACNA1G','CALY','CCDC81','CCKAR','CCL17','CCR8','CD72','CD79B','CENPA','CEP55','CETP','CLCA4','CLCN1','CLDN14','CLEC1A','CNGB3','CNKSR1','CNTFR','COQ3','CRH','CRX','CSHL1','CSRP3','CTSG','CXorf36','CYP2A7','DCC','DKKL1','DLX4','DNTT','DPEP3','DPYS','DSP','EDF1','EFNA2','ESPL1','FBRS','FBXO24','FCAR','FCN2','FGF8','FIP1L1','FLT3','FOXM1','FRS3','FSTL4','GABRA6','GCK','GDF10','GMIP','GNL2','GPR3','GPR4','GREB1','GRIK3','GTSE1','GYS2','H2AFX','HAMP','HCRTR2','HIST1H2BL','HIST1H2BM','HNRNPA0','HP1BP3','HPS4','HSPB6','HTR1B','HTR5A','IFNA1','IGLL1','IL12B','IMP4','IQCC','KCNJ13','KCNJ9','KCNQ4','KCNV2','KERA','KIF11','KIF14','KIF23','KIF4A','KLHL12','KNG1','KRI1','KRT12','KRT19','KRTAP1-3','LAMB4','LARP7','LILRB1','LILRB4','LLGL1','LRRTM4','LSM2','LTC4S','LY6H','MAG','MEN1','MEP1B','MKI67','MRPL12','MRPS15','MRTO4','MSTN','MUC6','MUSK','MYBL2','MYH4','MYH8','MYL7','NACA2','NEUROG1','NMUR1','NOL12','NOP56','NOTCH4','NUSAP1','NXF3','OMD','OMG','ONECUT1','OR7A5','OTOF','OTUD7B','P2RY14','PADI4','PAPOLB','PARK2','PCDH11Y','PCDHA5','PDE11A','PDE6D','PLA2G2D','PLK4','PMP2','PMPCA','POLA1','POU1F1','POU3F1','PPEF2','PPP4R2','PRB4','PRMT1','PRTN3','PSG11','PSORS1C2','PTTG1','QRSL1','R3HCC1','RAD23A','RAG2','RAPSN','RBP3','RFC2','RFC5','RNF25','RPE65','RPGRIP1','RRM2','SAC3D1','SCRT1','SGCA','SHCBP1','SIGLEC6','SIVA1','SLC12A1','SLC13A4','SLC25A31','SLC2A2','SLURP1','SMARCA4','SMC4','SNRPD1','SOX14','SPATA7','SPC25','SPTA1','SPTBN5','SRY','STIP1','STRN4','TACC3','TACR3','TCF3','TCL1B','TCOF1','TERT','TESK1','TGM3','THPO','TLX2','TNFRSF12A','TNNT2','TOP2B','TRA2A','TRPM1','TSSC1','TSSK2','TUT1','UBE2C','UCP3','VPREB1','VPREB3','VSX1','XPNPEP2','ZMYND10','ZNF155','ZNF407','ZNF674','ZNF81')
#DC
TGgene = c('ACHE','ADAM21','ADAMTS8','AHDC1','AKAP8L','ALCAM','ALDH1A2','ALOX15','ALOX15B','ANAPC2','ARL8B','ARPP21','ASPM','ATP1B4','BCL2L11','BCL2L13','BMP10','BTBD7','C1QA','C1QB','CA14','CALY','CAMK1G','CCDC81','CCL13','CCL17','CCL18','CCL19','CCL22','CCL23','CCL24','CCL8','CCR7','CD1A','CD1B','CD1C','CD1E','CD209','CD80','CD86','CD9','CEP350','CLCA4','CLDN14','CLEC10A','CLEC1A','CNGB3','CNKSR1','COQ3','CSRP3','CUL1','CXorf36','DKKL1','DNASE1L3','DPYS','EDF1','ETV3','F13A1','FBRS','FBXL4','FBXO24','FCER2','FGL2','FPR3','FSTL4','GMIP','GRIN1','GRSF1','GUCA1A','HAMP','HCRTR2','HIST1H2BL','HIST1H2BM','HK3','HLA-DQA1','HP1BP3','HPS5','HS3ST2','HSPB6','IL12B','IL21R','IMP4','IQCC','IRF4','KCNC3','KCNK13','KCNN1','KCNQ4','KCNV2','KERA','KIF23','KLHL12','KRI1','KRTAP1-3','LARP7','LILRB1','LILRB4','LOR','MAP3K13','MAP3K6','MCF2','MPHOSPH6','MRPL12','MRPS15','MS4A4A','MS4A6A','MYL7','NACA2','NAGPA','NECAP2','NFKB1','NMUR1','NOL12','NXF3','NXPH3','OR7A5','OTOF','PADI4','PAPOLB','PCDH11Y','PDE11A','PLD2','PLK4','PMPCA','PPP4R2','PRRG2','PSORS1C2','PTGES2','PTGIR','R3HCC1','RAB8A','RAPSN','RBP3','RNF2','RNF25','RPGRIP1','RRM2','RRP1B','SAMSN1','SIGLEC1','SLAMF1','SLAMF8','SLC13A4','SLC30A4','SLCO5A1','SLURP1','SNX11','SOX14','SPATA7','SPC25','SPINT2','SPTBN5','STAB1','STIP1','SUZ12','TACR3','TACSTD2','TBC1D13','TCF3','TDRD7','TESK1','TFEC','TGM3','THPO','TMEM131','TMSB10','TNFRSF12A','TNFRSF4','TNNT2','TRAF1','TREM2','TSSK2','TUT1','TXN','UBE2Z','UCP3','VAV2','VPREB1','VPREB3','VSX1','ZMYND10','ZNF155','ZNF407')
#iDC
TGgene = c('ALDH1A2','ALOX15','CCL13','CCL17','CCL18','CCL23','CCL24','CD209','CD86','CLEC10A','F13A1','FCER2','IL3RA','SPINT2')
#CD4+Tcell
TGgene = c('AAK1','ABCD2','ALG13','ANKRD55','APBB1','ARHGAP15','ARHGEF1','ARID5B','ASXL2','ATXN7','ATXN7L1','BAD','C14orf169','CA5B','CBLL1','CCNT2','CCR2','CCR4','CCR7','CCR8','CCR9','CD2','CD27','CD28','CD3D','CD3E','CD3G','CD4','CD40LG','CD5','CD6','CD7','CD96','CDC14A','CHMP7','CRLF3','CTLA4','CUBN','DDX31','DDX5','DIDO1','DLEC1','DNAH6','DNAJB1','DSC1','EZH1','FBXO11','FBXO21','FCN1','FNBP4','FOXP3','FUBP1','GIMAP4','GIMAP6','GOLGA4','GP5','GPR171','GPR183','HAUS3','HERC1','HIPK1','HMGN4','HMOX2','HNRNPU','ICOS','IL21R','IPCEF1','ITIH4','ITK','KLHL3','KRIT1','LAX1','LEPROTL1','LSM14A','LY9','MBNL1','MGAT2','MKL1','MORC2','MSL3','MTO1','NCK2','NFATC2IP','NFE2L2','NOL9','NR2C1','NUDCD3','OBSCN','OFD1','PABPC3','PHF3','PLCG1','PLCL1','PLXDC1','PNMA3','POLR3E','POU6F1','PPM1B','PPP2CA','PPWD1','PRMT2','PRPF38B','PSD','RAPGEF6','RBL2','RBM19','RBM25','RIC3','RSRC2','RXRG','SACM1L','SFPQ','SGSM2','SIRPG','SIT1','SLTM','SNPH','SNX19','SON','SPEG','SSTR3','SUPV3L1','THAP1','THUMPD1','TMEM123','TNFRSF4','TNFSF8','TOE1','TOMM20','TOR1AIP1','TPP2','TPT1','TRAF1','TRAF3IP3','TRAT1','TRIM46','TRMT2B','TSPYL1','TTN','TUBGCP5','TUG1','UBA3','UBASH3A','UBP1','UBQLN2','USP33','USP36','USP47','WBP11','ZAP70','ZBTB11','ZC3HAV1','ZCCHC11','ZFC3H1','ZNF335','ZNF611','ZNHIT6','ZXDC')
#CD4+naiveTcell
TGgene = c('09-Sep','AAK1','ACAP1','ACBD4','ANKRD55','APBB1','ARFRP1','ATXN7','BMS1','CA5B','CABIN1','CASP8','CCR7','CD2','CD226','CD247','CD27','CD28','CD3E','CD3G','CD4','CD40LG','CD5','CD6','CD7','CDK1','CDK10','CEPT1','CHMP7','CLC','CLUAP1','COPS7B','COQ6','CREBZF','CRLF3','CTLA4','CTSW','CUBN','CUL1','DDX31','DDX50','DIDO1','DNAJB1','DPEP2','DSC1','FAM193B','FCF1','FKTN','FNBP4','GIMAP6','GIN1','GLG1','GPSM3','GRAP2','GRK6','GZMM','HAUS3','HMOX2','ICOS','IDUA','IL16','INPP4A','INSL3','IPCEF1','ITK','JAK3','KDM3A','KLHL3','KRI1','KRT2','LAIR2','LEPROTL1','LIMD2','LY9','MAK','MKL1','MLH3','MLXIP','MSL3','MTRF1','NAA16','NCK2','NDFIP1','NOL9','NPAT','NUDCD3','NUDT9','NUMA1','NUP50','OBSCN','PACS1','PARP11','PHF1','PHF20L1','PIP4K2B','PLCG1','PLCL1','PLXDC1','POLR3E','POP5','PRMT2','PRMT3','RAB3GAP1','RAPGEF6','RBL2','RBMS1','REV1','RNF216','RNPEPL1','RPA3','RPAP2','RPL14','RPL38','RPLP2','RPRD2','RPS6','RXRB','SELPLG','SETD5','SH2B1','SIRPG','SIT1','SLTM','SNPH','SORCS3','STAP1','SUPV3L1','TATDN2','TBC1D5','TEX264','TMEM30B','TNFSF8','TNK1','TPP2','TRAF1','TRAF3IP3','TRAT1','TSPAN32','TUG1','UBASH3A','USP16','UTP20','VPS52','WDR6','WDR82','ZAP70','ZBTB40','ZNF263','ZNF264','ZNF394','ZNF609','ZNF76','ZNF780B')
#CD4+memoryTcell
TGgene = c('AAMP','ACD','ACTL6A','ADSL','AHCTF1','AKT2','AMBRA1','ANP32B','ANXA7','API5','ARHGAP15','ARL2','ARPC4','ATF1','ATG5','ATPIF1','ATXN10','AURKAIP1','BAG3','BCAS2','BTF3','BUB3','C11orf58','C12orf29','CBLL1','CBX3','CCNC','CCR4','CD2','CD226','CD28','CD2AP','CD3G','CD40LG','CD5','CD6','CD96','CDC40','CDK9','CDKN2AIP','CDV3','CEP57','CETN3','CLPX','CMPK1','CNBP','COPS4','COPS5','COX7C','CPSF6','CSNK1A1','CSNK2A2','CTBP1','CTLA4','CXCR6','DAD1','DAP3','DBF4','DDX3X','DDX50','DENR','DLEC1','DNAJA2','DNAJB1','DOHH','DPM1','DR1','EEF2','EID1','EIF2B5','EIF3E','EIF3L','EIF3M','EIF4G2','ERH','ESD','ETAA1','EXOC2','FARS2','FCF1','FNTA','FUBP3','FXR1','GABPA','GALR2','GATAD2A','GLOD4','GLUD1','GPR132','GPR15','GPR171','GPR183','GRPEL1','GZMA','GZMK','HDAC1','HINT1','HMGB2','HMGN4','HMOX2','HNRNPA0','HNRNPH3','HNRNPU','ICOS','IMP3','INTS8','ISCA1','ITK','JAK3','KARS','KBTBD4','KIF22','KTN1','LDHB','LIMS1','LIN7C','MAEA','MAGOH','MATR3','MEN1','METAP1','METTL5','MMADHC','MRPL11','MRPL20','MRPL44','MRPS18B','MRPS27','MRPS34','NAE1','NCL','NDUFS5','NKRF','NUPL2','OSGEP','PABPC4','PAPOLA','PCID2','PCNP','PDCD1','PDCD10','PFDN6','PFN1','PKP4','PLP2','POLD2','PPA2','PPID','PPIH','PPP1CB','PPP1CC','PPP2R5D','PPP6C','PREPL','PRPF18','PRPF19','PSMF1','PTGES3','PTPN11','PTPN4','RAD21','RANBP1','RANBP9','RBL2','RBM3','RBM34','RGS1','RNF34','RNF6','RPF1','RPL13','RPL13A','RPL36','RPL4','RPL5','RPL8','RPS19','RPS3','RPS6','RRP1B','RSL24D1','RUVBL1','RWDD1','SAFB2','SEC23IP','SERP1','SH2D1A','SLAMF1','SLC25A38','SLC25A6','SMAD2','SMC5','SMU1','SOD1','SP3','SPAG16','SRP9','SSNA1','STK16','SUB1','SUCLG1','SURF2','TBL3','THAP11','THOC7','THOP1','THRAP3','TINF2','TPP2','TPT1','TRA2A','TRAT1','TRMT112','TSN','TSSC1','TTC37','U2AF2','UBASH3A','UBE2D2','UBE2D3','UBE2N','UBIAD1','UBQLN2','UNC45A','UQCRC2','USP39','UXT','WDR46','ZC3H15','ZDHHC6','ZNF236','ZZZ3')
#CD4+Tcm
TGgene = c('09-Sep','AAK1','ACAP1','ADSL','AK1','ANAPC5','ANKRD55','APBB1','ARHGAP15','ARID5B','BAG3','BMPR1A','CA5B','CAMSAP1','CBLL1','CCR4','CCR6','CCR7','CCR8','CD2','CD247','CD28','CD4','CD40LG','CD48','CD5','CD6','CDC14A','CDKN2AIP','CNIH4','COL5A3','CORO7','CRY2','CSNK1D','CTLA4','DAB1','DDX24','DGCR14','DHX16','DIDO1','DLEC1','DNAI2','DNAJB1','DNMT1','DVL1','EDC4','EIF2B5','ERN1','EXOC1','FAM193B','FASTK','FBXL8','FXYD7','GGA3','GLTSCR2','GMEB2','GOLGB1','GP5','GPR15','GPR25','HNRNPUL1','HUWE1','ICOS','IDUA','IKZF1','IL2RA','INPP5E','IPCEF1','ITGB1BP1','ITIH4','ITK','JOSD1','KALRN','KBTBD2','KBTBD4','KDM3A','KRT1','LTA','LY9','MAN2C1','MCF2L2','MORC2','MSL3','MTO1','MYO16','NAP1L4','NCDN','NCK2','NDRG3','NFRKB','NME6','NPAT','NR2C1','NUP85','NXF1','OBSCN','OLAH','PCDHGA9','PIGG','PLCG1','PLCH2','PLCL1','PLXDC1','PNMA3','POLR2A','POU6F1','PSD','PSMD2','PURA','RAD9A','RAE1','RBM5','RPL38','RPS16','RPS21','RRS1','SARDH','SCNN1D','SGSM2','SIRPG','SLAMF1','SLC4A5','SMG5','SNPH','SNTG2','SOCS3','SORCS3','SPEG','SPTAN1','STK11','SUPT6H','TAB2','TCF20','TCF25','TFAP4','THAP11','TMEM30B','TNFRSF4','TNFSF8','TOMM7','TPO','TPR','TRADD','TRAF3IP3','TRAT1','TRIM46','TRMT61A','TSC1','UBASH3A','UBQLN2','USP10','USP36','USP39','USP4','WDR59','WDR6','XPC','YLPM1','ZC3HAV1','ZFYVE9','ZNF200','ZNF638','ZRSR2','ZXDB')
#CD4+Tem
TGgene = c('09-Sep','10-Sep','AHCTF1','AIRE','APBA3','ARAF','ARHGAP15','ASB6','BAG3','BCAS2','BRPF1','C19orf53','CABIN1','CAPZB','CCNT1','CCR2','CCR4','CCR6','CCR8','CCR9','CD2','CD28','CD3E','CD4','CD40LG','CD48','CD5','CD52','CD6','CDC14A','CDC37','CDC73','CGRRF1','CLIP1','CNOT1','COG4','COL5A3','COLQ','COPB1','COPS6','CTLA4','CXCR3','DCTN6','DDX56','DGCR14','DLEC1','DNAI2','DNAJB1','DYNC1H1','DYNLT1','E4F1','EIF3A','EIF3G','ELAC2','EMD','ERN1','ESYT1','FBXL8','FBXO31','FXYD7','FYCO1','GALNT8','GIPC1','GIT1','GLTSCR2','GOLGA7','GOLGB1','GORASP2','GPI','GPR15','GPR171','GPR25','GPSM3','GRM3','GZMK','HAUS3','HIST1H3A','HMGN4','IBTK','ICOS','IDE','IK','IL10RA','IL5RA','ITGB7','ITIH4','ITK','JMJD6','JTB','KBTBD4','KIAA0368','KLRB1','KRT1','LTA','LTK','LZTR1','MAML1','MAPKAPK5','MARK3','MCF2L2','MKNK1','MORC2','MTMR6','MUS81','MYO16','N4BP1','NCDN','NCK2','NDUFB2','NECAP2','NFKB1','NFRKB','NOL6','NSD1','OSM','PABPC3','PDCD1','PGK1','PLCL1','PSMD9','PVRIG','RABGGTA','RANBP9','RBL1','RBM10','RGS1','RIC8A','RIPK1','RNF7','RNPEPL1','RPL38','RPN2','RPS21','RRM1','RRS1','RWDD1','S100A11','SAMSN1','SEC24C','SELPLG','SF3A1','SF3B2','SIT1','SLAMF1','SLC35B1','SLC38A10','SLC9A6','SMAP1','SMARCC2','SOS1','SPSB3','SPTAN1','SRRT','SSR2','STUB1','STX17','TAF10','TAF1B','TBCC','TCF20','TCF25','THUMPD1','TKTL1','TNFRSF4','TNFSF8','TNIP2','TPO','TRADD','TRAF1','TRAPPC2L','TRAPPC4','TRAT1','TUFM','UBASH3A','UBIAD1','USP10','USP36','VPS54','XPC','ZC3HAV1','ZFYVE9','ZNF394')
#Treg
TGgene = c('ATG2B','BANP','CCR3','CCR4','CCR8','CD28','CD5','CTLA4','CXCR6','FOXP3','GALNT8','GPR25','HS3ST3B1','ICOS','IKZF4','IL10RA','IL2RA','IPCEF1','ITGB7','KCNA2','LAIR2','LAX1','LRP2BP','MCF2L2','MCM9','PLCL1','PPM1B','RGS1','SIT1','SPTAN1','STAM','TTN','TULP4','UBE4A','VPS54','ZCCHC8','ZFC3H1','ZMYM1','ZNF236')
#Th1
TGgene = c('CDC123','CHD1L','CHD4','COX10','CSTF1','CUEDC2','EIF2B2','FIBP','GNLY','HTRA2','IFNG','KIF20A','LAG3','MDC1','MNAT1','NCAPD3','NUP205','PKMYT1','POLD2','PPM1G','PSMD3','PTTG1','R3HDM1','RNPS1','RUVBL2','SLAMF1','SNRPC','TACO1','THOP1','TMEM39B','TRIM28','TTLL5','UBAP2','WDR18','WRAP53','ZBTB32')
#Th2
TGgene = c('BAG2','CDK2AP1','CEP55','CXCR6','GPR15','GZMA','GZMK','IL13','IL5','MAD2L1','NPHP4','NUP37','RAD50','RGS9','RNF34','RRAS2','RRM2','SLC25A44','SMAD2','THADA','TMEM39B','UBAP2')
#Plasma cells
TGgene = c('ABCB9','ACBD3','ACBD4','ADM2','AIPL1','ALG5','ALG9','AMPD1','APOA1','APOC3','ARHGEF16','ARL1','ARSA','ATF6','AUP1','AVIL','AVP','B4GALT3','BFSP2','BMP8B','C16orf58','C19orf73','C21orf2','C6orf25','CA7','CACNA1S','CAMP','CASP10','CAV1','CCDC121','CCDC33','CCDC40','CCDC88A','CCL25','CCNC','CD180','CD19','CD27','CD79A','CD79B','CDH15','CEACAM21','CELA2B','CHRNA4','CHRNG','CHST8','CLCNKB','CLINT1','CNKSR1','CNPY2','CNR1','CNTD2','COX6A2','CRYBA4','CRYBB1','CRYBB3','CRYGC','CSHL1','CSPP1','CUL7','CUX2','CYBA','CYP11A1','DAD1','DDN','DDOST','DEF8','DKKL1','DMTF1','DNAJC4','DNASE1L2','DOK3','DPAGT1','DRD4','DRD5','EBAG9','ELL','ENTPD1','EPO','ERGIC3','FBP2','FCRL2','FGF6','FKBP2','FN3K','FNDC3A','FTCD','GABRR2','GDF2','GH2','GLT8D1','GMPPA','GNB3','GOLGA3','GOLGA4','GOLGB1','GORASP2','GP9','GPLD1','GPR37L1','GPRC5D','GRIN1','GRM4','GRWD1','GUCA1A','HAND2','HDLBP','HEYL','HIST1H2BB','HOOK2','HSF4','HSP90B1','HSPA6','IDE','IFT52','IGF1','IMP4','IMP4AMPD1','IRF4','IRGC','ISCU','ITGA8','KCNN3','KCNQ4','KDELR2','KIAA0125','KLC2','KLF15','KNTC1','KRT10','LAX1','LBX1','LEFTY2','LMAN1L','LMAN2','LTB4R2','MAGEF1','MANF','MAPK8IP3','MARS','MAST1','MATN4','MBTPS1','MCTS1','MGAT2','MIA3','MIS12','MRPS31','MTDH','MYH13','MYL2','MYL7','NDOR1','NEUROG1','NGLY1','NOS2','NPAS1','NPPA','NPPC','NTRK1','OGFOD2','P2RY4','PABPC4','PCDHA5','PDE6A','PDIA2','PDIA6','PHOX2A','PICK1','PNOC','POMC','POU3F3','PPIB','PPIL2','PRDM14','PRDX4','PREB','PRM1','PRX','PTGER1','PTPRS','R3HCC1','RAB3A','RAD17','RALY','RASIP1','RAX','RFX2','RGS1','RGS13','RNF103','RNF113A','RNF208TNFRSF17','RNGTT','RPN1','RPN2','SAP30BP','SCFD1','SEC24A','SEC61A1','SEC61B','SEC61G','SEC62','SEC63','SEMA6C','SERP1','SHANK1','SHBG','SIPA1L3','SIX5','SLC13A2','SLC35B1','SLC35C2','SLC38A10','SLC5A2','SLC6A13','SLCO5A1','SMPD2','SNAPC4','SPATS2','SPCS1','SRP54','SRPRB','SS18','SSR1','SSR4','STMN4','SURF1','SYT5','T','TBL2','TBX4','TCF3','TERT','TG','THAP4','TIMM17B','TM9SF1','TMED10','TMED9','TMEM39A','TMEM59','TNFRSF13B','TNFRSF17','TNNT3','TP73','TRABD','TREH','TSHR','TSSK2','UBA5','UBE2G1','UBXN4','UFSP2','UGGT1','USP48','UTF1','VPREB1','VPREB3','VSX1','WDR45','WNT1','YIPF1','YIPF2','ZBP1','ZBP1ALPI','ZBP1AMPD1','ZDHHC4','ZNF133','ZNF142','ZNF37A','ZPBP')

for (i in TGgene){
  print(i)
  print(ensemblGenes[ensemblGenes$external_gene_name == i,'ensembl_gene_id'])
}

#pdf("PE24_TSNE_5000_naive_Tcell_external_gene_name_seurat.pdf")

exp =  BC_seurat@scale.data["ENSG00000058799",]
exp <- rbind(exp,exp)
exp[,] = 0

TGgene=c('ENSG00000150967','ENSG00000182827','ENSG00000181513','ENSG00000128165','ENSG00000120697',
  'ENSG00000086848','ENSG00000118137','ENSG00000130762','ENSG00000120805',
  'ENSG00000100299','ENSG00000118217','ENSG00000115307','ENSG00000135407','ENSG00000158850',
  'ENSG00000170819','ENSG00000140688','ENSG00000221916','ENSG00000160226',
  'ENSG00000003400','ENSG00000105974','ENSG00000176714',
  'ENSG00000141519','ENSG00000115355','ENSG00000112237','ENSG00000134061','ENSG00000177455',
  'ENSG00000139193','ENSG00000105369','ENSG00000007312','ENSG00000129910','ENSG00000007129',
  'ENSG00000113282','ENSG00000142675',
  'ENSG00000257727','ENSG00000118432','ENSG00000105219','ENSG00000100122',
  'ENSG00000104218','ENSG00000044090',
  'ENSG00000051523','ENSG00000140459','ENSG00000129562','ENSG00000244038','ENSG00000140995',
  'ENSG00000135164','ENSG00000110011','ENSG00000146094','ENSG00000172269',
  'ENSG00000069696','ENSG00000147654','ENSG00000105656','ENSG00000138185',
  'ENSG00000125991','ENSG00000132704','ENSG00000173486','ENSG00000167363',
  'ENSG00000102531','ENSG00000016864','ENSG00000144591',
  'ENSG00000111664','ENSG00000090615','ENSG00000144674','ENSG00000173230','ENSG00000115806',
  'ENSG00000112293','ENSG00000105447',
  'ENSG00000164107','ENSG00000115677','ENSG00000095066',
  'ENSG00000102878','ENSG00000166598','ENSG00000173110','ENSG00000119912','ENSG00000101052','ENSG00000017427',
  'ENSG00000136718','ENSG00000137265','ENSG00000136003','ENSG00000143603',
  'ENSG00000136240','ENSG00000174996','ENSG00000184445','ENSG00000186395',
  'ENSG00000122188','ENSG00000169223',
  'ENSG00000177383','ENSG00000145050','ENSG00000138834','ENSG00000166986','ENSG00000105613','ENSG00000124159',
  'ENSG00000140943','ENSG00000232119','ENSG00000168282','ENSG00000154305','ENSG00000167842','ENSG00000102738',
  'ENSG00000147649','ENSG00000188566',
  'ENSG00000151092','ENSG00000130751','ENSG00000163273',
  'ENSG00000111325','ENSG00000090621','ENSG00000132915','ENSG00000185615',
  'ENSG00000143870','ENSG00000100151','ENSG00000168081','ENSG00000115138',
  'ENSG00000166794','ENSG00000100023','ENSG00000123131','ENSG00000138073',
  'ENSG00000105227','ENSG00000105426','ENSG00000104679','ENSG00000105649','ENSG00000152942',
  'ENSG00000125970','ENSG00000105538','ENSG00000087903','ENSG00000090104','ENSG00000127074',
  'ENSG00000239305','ENSG00000125352','ENSG00000111880','ENSG00000163902','ENSG00000118705','ENSG00000161526',
  'ENSG00000092108','ENSG00000113615','ENSG00000058262','ENSG00000106803','ENSG00000132432','ENSG00000008952',
  'ENSG00000025796','ENSG00000120742','ENSG00000105738',
  'ENSG00000177045','ENSG00000121073','ENSG00000080189','ENSG00000157637','ENSG00000140675',
  'ENSG00000137571','ENSG00000135587','ENSG00000165684','ENSG00000123352','ENSG00000114902',
  'ENSG00000100883','ENSG00000144867','ENSG00000141380','ENSG00000124783','ENSG00000180879',
  'ENSG00000148290','ENSG00000129990','ENSG00000106638','ENSG00000071564',
  'ENSG00000176946','ENSG00000126768','ENSG00000100926','ENSG00000170348',
  'ENSG00000184840','ENSG00000176142','ENSG00000116209','ENSG00000240505','ENSG00000048462','ENSG00000130595',
  'ENSG00000078900','ENSG00000170638','ENSG00000165409','ENSG00000081307',
  'ENSG00000132388','ENSG00000144224','ENSG00000109775','ENSG00000136731','ENSG00000090686','ENSG00000171794',
  'ENSG00000128218','ENSG00000196998','ENSG00000058799',
  'ENSG00000130733','ENSG00000124256','ENSG00000136247','ENSG00000125846','ENSG00000115568','ENSG00000075407')

for (i in TGgene){
  test = BC_seurat@scale.data[i,]
  exp <- rbind(exp,test)
}

dim(exp)




exp <- exp[-c(1:2),]
save(exp, file="PE24_naiveCD8+T_ensg_cell_exp.RData")
load(file="PE24_naiveCD8+T_ensg_cell_exp.RData")

#PE24_sum_exp <- as.matrix(colSums(exp))

biocLite("vioplot")
library(vioplot)

pdf("PE24_DC_total_cell_boxplot.pdf")
#boxplot(PE24_sum_exp,notch=TRUE,col="darkgreen")
vioplot(PE24_sum_exp,PE25_sum_exp,PE26_sum_exp,PE29_sum_exp,names=c("PE24","PE25","PE26","PE29"),col=c("red","orange","darkgreen","blue"))
#stripchart(PE24_sum_exp,vertical=TRUE,method="jitter",add=TRUE,pch=1,col='red',cex=0.5)
dev.off()


TGgene='ENSG00000184271'
print(TGgene)
pdf(paste0("PE24_CD4+Tcell_",ensemblGenes[TGgene,"external_gene_name"],"_seurat.pdf"))
df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=BC_seurat@scale.data[TGgene,])
ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=1) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
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
  ) 
dev.off()



list=c('ENSG00000112208','ENSG00000111328','ENSG00000138180','ENSG00000172215','ENSG00000145649','ENSG00000113088','ENSG00000164109','ENSG00000131697','ENSG00000075188','ENSG00000113522','ENSG00000108370','ENSG00000170633','ENSG00000133818','ENSG00000171848','ENSG00000160785','ENSG00000175387','ENSG00000115970','ENSG00000121775','ENSG00000137073')
for (TGgene in list){
  print (TGgene)
  #pdf(paste0("test_",ensemblGenes[TGgene,"external_gene_name"],"_seurat.pdf"))
  df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=BC_seurat@scale.data[TGgene,])
  plot<-ggplot(df,aes(x=x, y=y, colour=expression)) + 
    geom_point(size=1) + 
    scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
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
    ) 
  ggsave(plot,file=paste0("PE24_PlasmaCell_",ensemblGenes[TGgene,"external_gene_name"],"_seurat.png"))
  #dev.off()
}


BC_seurat
colData(scesetFiltered)[,"total_counts"]

setwd("C:/Users/User/Dropbox/Angiogenesis_eunmin/aggr/outs/analysis/pca/10_components")
cellrangerPCA =read.csv("projection.csv")
eunmin_bacodes = sapply(colnames(scesetFiltered), function(x) substr(x, nchar(x)-15,nchar(x) ))
barcodes = sapply(cellrangerPCA["Barcode"], function(x) substr(x,1,16))
barcodes %in%  eunmin_bacodes
substr(cellrangerPCA["Barcode"], 1,16)

rownames(cellrangerPCA)

plot(cellrangerPCA[2:3])
df = data.frame(x=cellrangerPCA[2], 
                y=cellrangerPCA[3])
                # expression=exprs(scesetFiltered)[TGgene,])
ggplot(df,aes(x=x, y=y)) + 
  geom_point(size=1.0) 
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
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
  )

########################################### test cnv plot
PE24 = grepl(".*PE24.*", colnames(inferCNV_data))
table(PE24)
df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"],  
                expression=colSums(abs(inferCNV_data)))
ggplot(df,aes(x=x, y=y, colour=expression)) + 
  geom_point(size=1.0) + 
  scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
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
  ) 
dev.off()


df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=PE24)
centers = data.frame(df %>% dplyr::group_by(BC_seurat@ident) %>% summarize(x = median(x = x), 
                                                                           y = median(x = y)))
ggplot() + 
  geom_point(data=df, size=0.5, aes(x=x, y=y, colour=expression)) + 
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) +
  geom_text(data = centers,
            mapping = aes(x = x, y = y, label = BC_seurat.ident),
            size = 5)

setwd(data_path)

cancercell_list = read.csv("inferred_cancer_samples_2017_12_06_cancercells_5000.txt", sep="\t")
cancercell_list = apply(cancercell_list, 2, function(x) gsub("\\.","-",substr(x, 5, length(x))))
colnames_list = gsub("\\.","-",substr(colnames(inferCNV_data), 8, length(colnames(inferCNV_data))))
colnames(inferCNV_data)
inferCNV_data[cancercell_list %in% inferCNV_data]
colnames_list %in% cancercell_list

df = data.frame(x=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=BC_seurat@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=colnames_list %in% cancercell_list)
centers = data.frame(df %>% dplyr::group_by(BC_seurat@ident) %>% summarize(x = median(x = x), 
                                                                           y = median(x = y)))
ggplot() + 
  geom_point(data=df, size=0.5, aes(x=x, y=y, colour=expression)) + 
  ylab("Component 2") + 
  xlab("Component 1") + 
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(text = element_text(size=20),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), 
        axis.line=element_line(size=1),
        axis.ticks=element_line(size=1),
        legend.text=element_text(size=10), 
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x = element_text(size=20)
  ) +
  geom_text(data = centers,
            mapping = aes(x = x, y = y, label = BC_seurat.ident),
            size = 5)


#####################################################
#####infer CNV
#####################################################
# https://github.com/broadinstitute/inferCNV/wiki/ 참조
## R == 3.2.1 이어야만 함
# installation
#1 R 이용
library("devtools")
install_github("broadinstitute/inferCNV")
#2 source code download
R CMD install infercnv_0.1.tar.gz
# broadinstitute inferCNV download 

#input : exp matirx + genomic position file(optional)
## genomic position file 생성하는 python code
python gtf_to_position_file.py --attribute_name transcript_id your_reference.gtf your_gen_pos.txt
#convert_to_positional_file(args.input_gtf, args.output_positional, args.attribute_name)
# pip install [python_module]
##UCSC hg19 version whole gene (Homo_sapiens.GRCh37.75.gtf)
#Number of lines read: 2828312
#Number of comments: 5
#Number of entries: 63677
#Number of duplicate entries: 2764635
#Number of entries written: 63677

# 만들고 chr 순서대로 sorting하자 그래야 inferCNV x축이 순서대로 나옴
./scripts/inferCNV.R \
#  --ref_groups "1:235,236:280" \ #무시해도 되는 parameter
  --cutoff 4.5 \ 
  --noise_filter 0.3 \
  --output_dir quickstart \
##  --ref example/example_refs.txt \
  --vis_bound_threshold " -1,1" \
  example/example_expression.txt example/example_genomic_positions.txt


##### ESTIMATE : infer CNV 와 비슷한 tool
# http://bioinformatics.mdanderson.org/main/Tutorials
# installation
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

library(estimate)
#help(package="estimate")

setwd("/storage2/Project/CSC/RNA/PE24/03_quantification")

in.file <- file("PE24_RNA.genes.results","r")
out.file <- tempfile(pattern="estimate", fileext=".gct")
outputGCT(in.file, out.file)
filterCommonGenes(in.file, out.file)

