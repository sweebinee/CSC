setwd("/storage2/Project/CSC/RNA/03_Deconvolution/plot")

library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

toolname = 'DGIST_scRNA'
mtx <- read.table(file="DGIST_result_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)

mtx <- as.matrix(mtx)
mtx_melt <- melt(mtx[,3:8])
colnames(mtx_melt) <- c('cellTypes','sample','value')

ggheatmap = ggplot(mtx_melt, aes(sample, cellTypes, fill = value)) + 
labs(title = toolname)+
geom_tile(color = "black") + 
geom_text(aes(label=round(value,2)), color='black') +
scale_fill_gradient2(low = "white", high = "red", limit = c(0, 1), space = "Lab", name = "cell type freq") + 
theme_minimal(base_size = 12, base_family = "") + 
theme(axis.text.x=element_text(angle=45,vjust=1,size=15,hjust=1),axis.text.y=element_text(vjust=1,size=15,hjust=1))+

png(paste0(toolname,"_CRITERIA_01_heatmap.png"),width=600, height=600)
ggheatmap
dev.off()

###############################################################################
#scatter plot
library(plotly)
library(ggpubr)


mtx <- read.table(file="DGIST_result_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)
mtx <- as.matrix(mtx)
mtx_melt <- melt(mtx[,3:8])
#> mtx_melt
#                        Var1 Var2        value
#1                     B_cell PE24 0.0610323672
colnames(mtx_melt)[3] <- "scRNA"

toolnames = c('MCPcounter','CIBERSORT','Xcell')
for(i in toolnames){
  assign(i, read.delim(paste0("01CRITERIA_",i,".txt"), header=T, stringsAsFactors=F, row.names=1))
  assign(paste0(i,"_melt"),melt(as.matrix(get(i))))
  assign(paste0(i,"_merge"),merge(mtx_melt,get(paste0(i,"_melt")),by=c('Var1','Var2')))
  png(paste0(i,"_CRITERIA_01_dot.png"),width=600, height=600)
  plot(x=get(paste0(i,"_merge"))$scRNA,y=get(paste0(i,"_merge"))$value)
  dev.off()
}
#only CIBERSORT
CIBERSORT <- read.delim("mark05_CIBERSORT.txt", header=T, stringsAsFactors=F, row.names=1)
CIBERSORT <- melt(as.matrix(CIBERSORT))
mtx_melt <- merge(mtx_melt,CIBERSORT,by=c('Var1','Var2'))
colnames(mtx_melt)[4]<-"CIBERSORT"
merge_mtx <- mtx_melt


colnames(MCPcounter_merge)[4] <- "MCPcounter"
merge_mtx <- merge(MCPcounter_merge,CIBERSORT_melt,by=c('Var1','Var2'))
colnames(merge_mtx)[5]<-"CIBERSORT"
merge_mtx <- merge(merge_mtx,Xcell_melt,by=c('Var1','Var2'))
colnames(merge_mtx)[6]<-"Xcell"


# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y,method="spearman"), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
png("mark04_dot.png",width=600, height=600)
pairs(merge_mtx[,3:4], 
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      main="scRNA vs CIBERSORT\n")
dev.off()


#scRNA-CIBERSORT, scRNA-Xcell only scatter plot
library(ggrepel)

for(i in 1:nrow(merge_mtx)){
  rownames(merge_mtx)[i]<-paste0(merge_mtx$Var1[i],":",merge_mtx$Var2[i])
}

#remove Tcell
merge_mtx<-merge_mtx[-c(79:84),]

# Add text to the plot
.labs <- rownames(merge_mtx)
b <- ggscatter(merge_mtx, x = "scRNA", y = "CIBERSORT",
  add = "reg.line", 
#   add.params = list(color = "blue", fill = "lightgray"),
  color = "Var1", 
  palette=c("#FF0000", "#FF5500","#FFAA00","#FFFF00","#AAFF00","#00FF2B","#00FFD4","#00D4FF","#00AAFF","#0055FF","#5500FF","#AA00FF","#FF00AA","#FF0055"),
  conf.int = TRUE)

png("mark05_test.png",width=800, height=600)
b + xlim(0.001, 1)+ylim(0.001, 1)+
  geom_point(aes(color = Var1)) +
  geom_smooth(method='lm',se = FALSE, fullrange = TRUE)+
  stat_cor(label.x = 0.003,method="spearman")+
  geom_text_repel(aes(label = .labs,  color = Var1), size = 3)+
  scale_color_manual(values = c("#FF0000", "#FF5500","#FFAA00","#FFFF00","#AAFF00","#00FF2B","#00FFD4","#00D4FF","#00AAFF","#0055FF","#5500FF","#AA00FF","#FF00AA","#FF0055"))+
  labs(title = "scRNA vs CIBERSORT\n", color = "cellTypes\n")
dev.off()

#spearman correlation
cor(merge_mtx$scRNA, merge_mtx$CIBERSORT, method="pearson")
cor(merge_mtx$CIBERSORT, merge_mtx$scRNA, method="spearman")

######################################################################
#version 2 in the local
######################################################################
library(plotly)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(ggthemes)

main_dir = "/home/subin/Desktop/CSC"
setwd(main_dir)

file_dir = paste0(main_dir,"/25Mar19/")
ver = "mark06"

mtx <- read.table(file="DGIST_result_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)
mtx <- as.matrix(mtx)
mtx_melt <- melt(mtx[,3:8])
colnames(mtx_melt)[3] <- "scRNA"

CIBERSORT <- read.delim(paste0(ver,"_CIBERSORT.txt"), header=T, stringsAsFactors=F, row.names=1)
CIBERSORT <- melt(as.matrix(CIBERSORT))
merge_mtx <- merge(mtx_melt,CIBERSORT,by=c('Var1','Var2'))
colnames(merge_mtx)[4]<-"CIBERSORT"

#scRNA-CIBERSORT, scRNA-Xcell only scatter plot
for(i in 1:nrow(merge_mtx)){
  rownames(merge_mtx)[i]<-paste0(merge_mtx$Var1[i],":",merge_mtx$Var2[i])
}

###################
for(CT in unique(merge_mtx$Var1)){
plot_m<-merge_mtx[merge_mtx$Var1==CT,]

ggscatter(plot_m, color="black", x = "scRNA", y = "sumSC",xlab=FALSE, ylab=FALSE,size=1,
  ggtheme=theme_classic(base_size=10)) +
  coord_fixed()+ #x, y axis 1:1
  coord_cartesian(expand = FALSE,xlim=c(0,1),ylim=c(0,1))+ # 여백제거
  geom_abline(intercept=0, slope=1, color="grey")+
  annotate(geom="text", x=0.2, y=0.9, label=paste0("R=",round(cor(plot_m$sumSC, plot_m$scRNA, method="pearson"),2)), color="black",fontface=2)
ggsave(paste0(file_dir,"cluster_mark03_scRNA_sum.png"),width = 2, height = 2)
}

###################
plot_m<-merge_mtx[merge_mtx$scRNA>0&merge_mtx$CIBERSORT>0,]
b <- ggscatter(plot_m, x = "scRNA", y = "CIBERSORT",
  add = "reg.line", 
  #add.params = list(color = "blue", fill = "lightgray"),
  color = "Var1" 
  #  conf.int = TRUE #주변에 영역표시
)
# Add text to the plot
.labs <- rownames(plot_m)
color_pl <- c("B_cell"="#FF0000", "Bone_marrow_cells"="#FF5500","CMP"="#FFAA00","Dendritic"="#FFFF00",
  "Erythroblast"="#AAFF00","GMP"="#00FF2B","Haematopoietic_stem_cells"="#00FFD4","Macrophage"="#00D4FF",
  "Monocyte"="#00AAFF","Myelocyte"="#0055FF","Neutrophil"="#5500FF","NK_cell"="#AA00FF",
  "Pro-Myelocyte"="#FF00AA","T_cell"="#FF0055")

pdf(paste0(file_dir,ver,"_cor.pdf"),width=12,height=10)
b +xlim(0, 1)+ylim(0, 1)+
  geom_point(aes(color = Var1)) + 
  theme_hc(style = "darkunica") + 
  theme(legend.position="right") + 
  #scale_colour_hc("darkunica") +
  #geom_smooth(method='lm',se = FALSE, fullrange = TRUE)+ #전체 cor line
  geom_text_repel(aes(label = .labs,  color = Var1), size = 3)+
  scale_color_manual(name="cellTypes",values = color_pl)+
  labs(title = paste0("scRNA vs CIBERSORT\n",ver), color = "cellTypes\n")+
  #total cor
  annotate(geom="text", x=0.7, y=0.97, label=paste0("R=",round(cor(merge_mtx$CIBERSORT, merge_mtx$scRNA, method="pearson"),2)), color="white",fontface=2)+
  annotate(geom="text", x=0.7, y=0.93, label=paste0("B_cell : R= ",round(plot_cor["B_cell",],2)), color=color_pl[1])+
  #annotate(geom="text", x=0.7, y=0.91, label=paste0("Bone_marrow_cells : R= ",round(plot_cor["Bone_marrow_cells",],2)), color=color_pl[2])+
  #annotate(geom="text", x=0.7, y=0.89, label=paste0("CMP : R= ",round(plot_cor["CMP",],2)), color=color_pl[3])+
  annotate(geom="text", x=0.7, y=0.87, label=paste0("Dendritic : R= ",round(plot_cor["Dendritic",],2)), color=color_pl[4])+
  #annotate(geom="text", x=0.7, y=0.85, label=paste0("Erythroblast : R= ",round(plot_cor["Erythroblast",],2)), color=color_pl[5])+
  annotate(geom="text", x=0.7, y=0.83, label=paste0("GMP : R= ",round(plot_cor["GMP",],2)), color=color_pl[6])+
  annotate(geom="text", x=0.7, y=0.81, label=paste0("Haematopoietic_stem_cells : R= ",round(plot_cor["Haematopoietic_stem_cells",],2)), color=color_pl[7],fontface=2)+
  annotate(geom="text", x=0.7, y=0.79, label=paste0("Macrophage : R= ",round(plot_cor["Macrophage",],2)), color=color_pl[8],fontface=2)+
  annotate(geom="text", x=0.7, y=0.77, label=paste0("Monocyte : R= ",round(plot_cor["Monocyte",],2)), color=color_pl[9])+
  #annotate(geom="text", x=0.7, y=0.75, label=paste0("Myelocyte : R= ",round(plot_cor["Myelocyte",],2)), color=color_pl[10])+
  #annotate(geom="text", x=0.7, y=0.73, label=paste0("Neutrophil : R= ",round(plot_cor["Neutrophil",],2)), color=color_pl[11])+
  annotate(geom="text", x=0.7, y=0.71, label=paste0("NK_cell : R= ",round(plot_cor["NK_cell",],2)), color=color_pl[12])+
  #annotate(geom="text", x=0.7, y=0.69, label=paste0("Pro-Myelocyte : R= ",round(plot_cor["Pro-Myelocyte",],2)), color=color_pl[13])+
  annotate(geom="text", x=0.7, y=0.67, label=paste0("T_cell : R= ",round(plot_cor["T_cell",],2)), color=color_pl[14],fontface=2)
dev.off()

#spearman correlation
plot_cor <- as.data.frame(matrix(,nrow=14,ncol=1))
rownames(plot_cor) <- unique(merge_mtx$Var1)
colnames(plot_cor) <- "corr"
for(CT in unique(merge_mtx$Var1)){
  c<-cor(merge_mtx[merge_mtx$Var1==CT,"scRNA"], merge_mtx[merge_mtx$Var1==CT,"CIBERSORT"], method="pearson")
  plot_cor[CT,1]<-c
}
cor(merge_mtx$CIBERSORT, merge_mtx$scRNA, method="spearman")

######################################################################
#version 3 in the local : cluster ver.
######################################################################
library(plotly)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(ggthemes)

main_dir = "/home/subin/Desktop/CSC"
#MAC OS 
#main_dir = "/Users/subin/Desktop/CSC"

setwd(main_dir)

file_dir = paste0(main_dir,"/25Mar19/")
ver = "cluster_mark03"

mtx <- read.table(file="DGIST_cluster_percentage.txt",sep='\t',header=TRUE,stringsAsFactor=FALSE,row.names=1)
mtx_melt <- melt(as.matrix(mtx))
#mtx_melt <- melt(mtx[,3:8])
colnames(mtx_melt)[3] <- "scRNA"

sumSC <- read.delim(paste0("scRNA_SUM_CIBERSORT.txt"), header=T, stringsAsFactors=F, row.names=1)
sumSC <- melt(as.matrix(sumSC))
merge_mtx <- merge(mtx_melt,sumSC,by=c('Var1','Var2'))
colnames(merge_mtx)[4]<-"sumSC"

CIBERSORT <- read.delim(paste0(ver,"_CIBERSORT.txt"), header=T, stringsAsFactors=F, row.names=1)
CIBERSORT <- melt(as.matrix(CIBERSORT))
merge_mtx <- merge(merge_mtx,CIBERSORT,by=c('Var1','Var2'),all.x=TRUE)
colnames(merge_mtx)[5]<-"CIBERSORT"
merge_mtx<-merge_mtx[order(merge_mtx[,"Var1"]),]

merge_mtx$Var1 <- as.character(merge_mtx$Var1)

#scRNA-CIBERSORT, scRNA-Xcell only scatter plot
for(i in 1:nrow(merge_mtx)){
  rownames(merge_mtx)[i]<-paste0(merge_mtx$Var1[i],":",merge_mtx$Var2[i])
}

plot_m<-merge_mtx[merge_mtx$scRNA>0&merge_mtx$CIBERSORT>0,]

plot_m <- plot_m[plot_m$Var1%in%c("1","4","8","9","10","14","15","16","21"),]

b <- ggscatter(plot_m, x = "scRNA", y = "CIBERSORT",
  add = "reg.line", 
  color = "Var1" )

.labs <- rownames(plot_m)

#spearman correlation
plot_cor <- as.data.frame(matrix(,nrow=22,ncol=1))
rownames(plot_cor) <- unique(merge_mtx$Var1)
colnames(plot_cor) <- "corr"
for(CT in unique(merge_mtx$Var1)){
  c<-cor(merge_mtx[merge_mtx$Var1==CT,"scRNA"], merge_mtx[merge_mtx$Var1==CT,"CIBERSORT"], method="pearson")
  plot_cor[CT,1]<-c
}

pdf(paste0(file_dir,ver,"_cor.pdf"),width=12,height=10)
b +xlim(0, 1)+ylim(0, 1)+
  geom_point(aes(color = Var1)) + 
  theme_hc(style = "darkunica") + 
  theme(legend.position="right") + 
  #scale_colour_hc("darkunica") +
  #geom_smooth(method='lm',se = FALSE, fullrange = TRUE)+ #전체 cor line
  geom_text_repel(aes(label = .labs,  color = Var1), size = 3)+
  labs(title = paste0("scRNA vs CIBERSORT\n",ver), color = "cellTypes\n")
dev.off()

##########################3
#bar plot으로 그려보자
library(ggplot2)
library(scales)
library(gridExtra)


sample=c("PE23_1","PE23_2","PE24","PE25","PE26","PE29","PE32","PE36")
cluster=unique(merge_mtx$Var1)
per_mtx <- as.data.frame(matrix(, nrow=264, ncol=4))
colnames(per_mtx) <- c("PE","CLUSTER","TOOL","PERCENTAGE")
i = 1
for(PE in sample){
  for(CL in cluster){
    for(tool in c("scRNA","sumSC","CIBERSORT")){
      per_mtx[i,"PE"] <- PE
      per_mtx[i,"CLUSTER"] <- CL
      per_mtx[i,"TOOL"] <- tool
      rn <- paste0(CL,":",PE)
      per_mtx[i,"PERCENTAGE"] <- merge_mtx[rn,tool]
      i = i+1
    }
  }
}

fill <- c("#0073C2FF","#99CC00","#EFC000FF")

par(mfrow = c(6, 1))

for(i in sample){
  dt<-per_mtx[per_mtx$PE==i,]
    ggplot(data=dt,aes(x=factor(as.numeric(dt$CLUSTER)), 
      y=dt$PERCENTAGE, fill=factor(dt$TOOL)))+
    coord_cartesian(ylim=c(0,1))+
    labs(x="", y="") +
    geom_bar(stat="identity", position = position_dodge(),width = 0.7)+
    #ggtitle(paste0(i,"\n Composition of clusters (%)")) +
    #theme_minimal() +
    scale_fill_manual(values=fill) +
    #guides(fill=FALSE)+
    #theme(legend.position="bottom", legend.direction="horizontal",
    #      legend.title = element_blank()) +
    theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
    theme(axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10))
  ggsave(paste0(file_dir,"cluster_mark03_plot_",i,".png"),width = 20, height = 5)
}


##
#what is different bw bulk RNA-seq gene vs scRNA-seq clusters marker gene?
bulk <- read.delim("/storage2/Project/CSC/RNA/03_Deconvolution/CSC_RNA_TPM_rm0_HUGO_rmdup.txt",sep='\t',header=TRUE,row.names=1,stringsAsFactors=FALSE)
bulk_gene <- rownames(bulk)
> length(bulk_gene)
[1] 28221
> length(union_DEG)
[1] 1100
> length(intersect(union_DEG,bulk_gene))
[1] 1087
