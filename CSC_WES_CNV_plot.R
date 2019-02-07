setwd("/storage2/Project/CSC/WES/04_CNV/EXCAVATOR2/w_10K")

# loading annotated CNV files
sample=c('PE24','PE32','PE36','PE25','PE26','PE29')

for(i in sample){
  assign(i, read.delim(paste0(i,"_EX_w10K.txt"), header=T, stringsAsFactors = F))
}

gene=NULL
for(i in sample){
  samplegene=get(i)$HGNC.symbol
  gene=unique(c(gene, samplegene))
}
gene<-sort(gene)

#  ------------------------------------------------------------------------

## make binary matrix
cnvmat <- matrix(NA, nrow = length(gene), ncol=length(sample))
colnames(cnvmat)=sample
rownames(cnvmat)=gene
cnvmat <- cnvmat[-1,]

for(i in sample){
  for(j in rownames(cnvmat)){
    samplemut <- unique(get(i)$HGNC.symbol)
    if(j %in% samplemut)
      cnvmat[j,i]=1
    else
      cnvmat[j,i]=0
  }
}

> dim(cnvmat)
[1] 1434    6

## compare with DGIST BCSC-like cluster marker 
BCSCmarker = c('MTATP6P1','MT-ND4L','MT-TP','MT-ND6','XIST','KCNQ1OT1','MTRNR2L1','MT-ND5','SPTBN1','AC027290.2','EIF3B','TRIM28','PRRC2A','RCC2','FLNB','MAP7D1','PTPRF','HUWE1','KHDC4','FAM120A','WEE1','FLNA','MT-ND3','FASN','AKAP13','RRBP1','MT-ND1','TFRC','DUS1L','LYPD3','BCYRN1','TSC22D2','AC023157.3','NEAT1','SLC7A5','CHERP','TAPBP','AFF4','IMPAD1','MT-RNR1','ALYREF','ATP2A2','MT-CYB','FAM84B','CAND1','COLGALT1','CDH1','MT-ND2','ARID1A','MT-ATP6','PABPN1','UBE2H','TENT5A','NKTR','KDM2A','NPTN','KIF1B','KMT2C','MT-CO1','MT-ND4','MT-CO2','NCKAP1','SLC38A2','MT-RNR2','SPTAN1','MAT2A','MT-CO3','DSG2','RHOB','CD81','IRF2BP2','MYH9','RYBP','TRPS1','CLIP1','RABL6','H1FX','PPP1R14B','MKNK2','DYNC1H1','ENAH','TCF7L1','EIF4G1','KMT2E','NFIC','CLTC','CTSZ','AKAP9','RAP2B','PITRM1','IRF1','CEBPD','ADAM15','BRD2','PTPRK','VEGFA','TMEM165','CCNL1','ANKRD54','ATP1B1','UBE2S','WSB1','PIM3','SF1','C6orf132','SON','PNISR','PPP1R15A','SUPT5H','IFRD1','GPRC5A','SIRT7','PRRC2C','COL6A1','PRKDC','ZFP36L1','SYNGR2','SLC52A2','SRM','HDGF','EDN1','VMP1','KPNB1','GNB1','AL118516.1','LSR','GADD45B','PRPF38B','C3','DSP','CYR61','HK2','TRA2B','CLDN4','MALAT1','CTSB','GIPC1','CAPNS1','RBM25','JUP','NOC2L','SOX4','SYNE2','SDC4','RBM39','SLC9A3R1','MLF2','EIF3A','XBP1','TMED9','EFNA1','AUP1','ACTN4','MTDH','KLF6','SDC1','ATP1A1','NDUFB9','HSPA5','CALM2','GAL','FUS','GNAS','HES1','WARS','GBP1','STMN1','GADD45A','JUN','IFI6','RPS6','YWHAZ','PEG10','HNRNPH1','ELF3','COX6A1','CALD1','RPL31','RPS3A','EEF1A1','MFGE8','RPLP0','GAPDH','TUBB','PTMA')
final <- cnvmat[rownames(cnvmat)%in%BCSCmarker,]

#> dim(final)
#[1] 351   6


#  ------------------------------------------------------------------------
## draw heatmap

library(gplots)

mycolor <- colorRampPalette(c("white", "black"))(n = 10)

png("CSC_WES_CNV_BCSCmarker_plot.png",width=800,height=800)
heatmap.2(final, Rowv=F, Colv=F, dendrogram="none", 
          col=mycolor,trace="none", key=F,
          xlab = "sample",
          ylab ="gene",
          margin=c(10,15), 
		#cexRow=1, cexCol=1, cex.main=0.5
          main="BCSC-like cluster marker CNV")
dev.off()

