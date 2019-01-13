setwd("/storage2/Project/CSC/RNA/03_Deconvolution")

cut -f1,5,6,7 /storage2/Project/CSC/RNA/02_Preprocessing/03.RSEM/PE36.genes.results > /storage2/Project/CSC/RNA/03_Deconvolution/PE36.txt
#CSC_RNA_rawCount.txt
#CSC_RNA_TPM.txt
#CSC_RNA_FPKM.txt
sed '/\t0.00\t0.00\t0.00\t0.00/d' CSC_RNA_rawCount.txt > CSC_RNA_rawCount_rm0.txt

load("/storage2/Project/CSC/10X/DGIST_data/ensemblGenes2018-11-16.RData")
df <- read.table("CSC_RNA_rawCount_rm0.txt", sep="\t", header=TRUE)
genes <- df$gene_id
#G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart=ensemblGenes)
df_merged <- merge(df,ensemblGenes,by.x="gene_id",by.y="ensembl_gene_id")
df_exist_gene_symbol <- df_merged[df_merged$external_gene_name!="",]
df_final <- df_exist_gene_symbol[,c('external_gene_name','PE24','PE25','PE26','PE29','PE32','PE36')]

#write.table(df_exist_gene_symbol, file = "CSC_RNA_rawCount_rm0_HUGO_raw.txt", sep="\t", row.names=FALSE)
#write.table(df_final, file = "CSC_RNA_rawCount_rm0_HUGO.txt", sep="\t", row.names=FALSE)
library(plyr)
dup_gene <- subset(count(df_final,vars=c("external_gene_name")),freq>1)$external_gene_name
rmdup <- df_final[!(df_final$external_gene_name%in%dup_gene),]
for(i in dup_gene){
  dup <- df_final[df_final$external_gene_name==i,]
  test <- as.data.frame(rbind(c(i,apply(dup[,-1],2,mean))))
  colnames(test)[1] <- "external_gene_name"
  rmdup <- rbind(rmdup,test)
}
#write.table(rmdup, file = "CSC_RNA_rawCount_rm0_HUGO_rmdup.txt", sep="\t", row.names=FALSE)

#####################################################################
###MCPcounter
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("")
#"devtools","curl","ggplot","reshape"

library(devtools)
install_github("ebecht/MCPcounter",ref="master", subdir="Source")

################################################################################################
MCPcounter.estimate         package:MCPcounter         R Documentation

MCPcounter.estimate

Description:
     this function produces a matrix with abundance estimates from an
     expression matrix

Usage:
     MCPcounter.estimate(expression,featuresType=c("affy133P2_probesets","HUGO_symbols","ENTREZ_ID")[1],
                     probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Si
gnatures/probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                     genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signat
ures/genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
     )
     
Arguments:

expression: matrix or data.frame with features in rows and samples in
          columns
featuresType: type of identifiers for expression features. Defaults to
          "affy133P2_probesets" for Affymetrix Human Genome 133 Plus
          2.0 probesets. Other options are "HUGO_symbols" (Official
          gene symbols) or "ENTREZ_ID" (Entrez Gene ID)

probesets: Data frame of probesets transcriptomic markers and
          corresponding cell populations. Fetched from github by a call
          to read.table by default, but can also be a data.frame

   genes: Data frame of genes transcriptomic markers (HUGO symbols or
          ENTREZ_ID) and corresponding cell populations. Fetched from
          github by a call to read.table by default, but can also be a
          data.frame

Value:

     matrix with cell populations in rows and samples in columns

Note:

     This is a contribution from the Tumor Identity Cards (CIT) program
     founded by the 'Ligue Nationale Contre le Cancer' (France): <URL:
     http://cit.ligue-cancer.net>. For any question please contact
     <URL: CITR@ligue-cancer.net>

Author(s):

     Etienne Becht

Examples:

     ExampleEstimates=MCPcounter.estimate(MCPcounterExampleData,featuresType="affy133P2_probesets")
     heatmap(as.matrix(ExampleEstimates),col=colorRampPalette(c("blue","white","red"))(100)) 
     
############################################################################################################


cat CSC_RNA_rawCount_rm0_HUGO.txt | cut -f1 | sort | uniq -c | sort
cat CSC_RNA_rawCount_rm0_HUGO.txt | grep -E 'U1|U3|Metazoa_SRP' > pe24_duplicate.txt
sed '/7SK\t0/d' CSC_rawCount_rm0_HUGO.txt|sed '/ALG1L9P\t/d'|sed '/BMS1P4\t/d'|sed '/COG8\t/d'|sed '/CYB561D2\t/d'|sed '/DGCR5\t/d'|sed '/DNAJC9-AS1\t/d'|sed '/GOLGA8M\t/d'|sed '/H2BFS\t/d'|sed '/LINC01238\t/d'|sed '/LINC01297\t/d'|sed '/LINC01422\t/d'|sed '/LINC01481\t/d'|sed '/MATR3\t/d'|sed '/POLR2J4\t/d'|sed '/PRSS50\t/d'|sed '/RABGEF1\t/d'|sed '/RGS5\t/d'|sed '/RNA5-8S4\t/d'|sed '/SNORA12\t/d'|sed '/SNORA63\t/d'|sed '/SNORD22\t/d'|sed '/SPATA13\t/d'|sed '/STPG4\t/d'|sed '/TMSB15B\t/d'|sed '/Y_RNA\t/d'|sed '/uc_338\t/d'|sed '/U3\t0/d'|sed '/Metazoa_SRP\t/d'|sed '/U3\t1\t2\t1\t2/d'|sed '/U1\t0\t2\t0\t1/d'|sed '/U1\t2\t1\t1\t0/d'|sed '/U1\t0\t0\t0\t1/d'|sed '/U1\t0\t4\t5\t0/d'|sed '/EMG1\t0\t24.21/d'|sed '/EMG1\t159.63/d'|sed '/SCO2\t80.8/d'|sed '/SCO2\t18.72/d' > 1.txt
cat 1.txt duplicate_mean.txt > CSC_rawCount_rm0_rmdup_HUGO.txt

setwd("/storage2/Project/CSC/RNA/04_MCPcounter")
setwd("/home/subin95/CSC/MCPcounter")

library(MCPcounter)

exp <- read.table("/home/subin95/CSC/RNA/CSC_rawCount_rm0_rmdup_HUGO.txt", sep="\t", header=TRUE,row.names=1)

CSC=MCPcounter.estimate(exp,featuresType="HUGO_symbols")
heatmap(as.matrix(CSC),sideColors=c("darkgreen", "yellowgreen"),col=colorRampPalette(c("blue","white","red"))(100)) 


library(ggplot2)
library(reshape)
library(plyr)
library(scales)
library(TDMR)

exp <- read.table("/home/subin95/CSC/MCPcounter/MCPresult.txt", sep="\t", header=TRUE)
exp.m <- melt(exp)
exp.m <- ddply(exp.m, .(variable), transform,rescale = rescale(value))
colnames(exp.m)[1]<-"Cell_Type"
colnames(exp.m)[2]<-"Samples"


pdf("MCP.pdf")
(p <- ggplot(exp.m, aes(Samples, Cell_Type)) 
	+ geom_tile(aes(fill = rescale),colour = "white") 
	+ scale_fill_gradient2(low = "blue",mid = "white",high = "red")
	+ geom_text(aes(label = round(value, 1))) )
base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "Samples", y = "Cell Types") + scale_x_discrete(expand = c(0, 0)) +scale_y_discrete(expand = c(0, 0)) 
#+ opts(legend.position = "none",axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

dev.off()
