
setwd('/storage2/Project/CSC/10X/DGIST_data')

DGIST=readRDS(file="BC_seurat.rds")

DGIST@raw.data      DGIST@var.genes     DGIST@meta.data     DGIST@assay         DGIST@cell.names    DGIST@calc.params   DGIST@misc
DGIST@data          DGIST@is.expr       DGIST@project.name  DGIST@hvg.info      DGIST@cluster.tree  DGIST@kmeans        DGIST@version
DGIST@scale.data    DGIST@ident         DGIST@dr            DGIST@imputed       DGIST@snn           DGIST@spatial    

#####################################################

library(Seurat)
library(scater)
library(dplyr)

pdf("DGIST_tsne.pdf", width=15, height=9)
TSNEPlot(DGIST, pt.size = 0.5, label.size = 4, do.label = TRUE, cex=2)
dev.off()

cell.num <- table(DGIST@ident)


library(ggplot2)
library(reshape)
library(plyr)
library(scales)
library(TDMR)

setwd("/home/subin95/CSC/10X")
df <- read.table("DGIST_cell.txt", sep="\t", header=TRUE)
df.m <- melt(df)

df_rmSum <- df.m[df.m$variable!="sum",]
df_final <- cbind(df_rmSum, percentage=c("1.73%","30.97%","6.07%","7.27%","53.96%","1.49%","39.11%","2.95%","9.65%","46.76%","7.4%","8.95%","3%","1.46%","76.59%","10%","7.24%","7.75%","0.94%","73.74%","2.06%","40.03%","8.26%","10.09%","39.54%"))
write.table(df_final, file = "DGIST_cell.txt", sep="\t", row.names=FALSE)

df <- read.table("DGIST_cell.txt", sep="\t", header=TRUE)

pdf("DGIST_cell.pdf", width=10, height=10)
p <- ggplot(df, aes(fill=Cell.Type, y=value, x=variable))+geom_bar( stat="identity")+geom_text(aes(label = percentage),size =3.5, position = position_stack(vjust = 0.5),colour="white",fontface = "bold")
base_size <- 15
p + ggtitle("DGIST # of Cells by Cell Types") + 
theme_grey(base_size = base_size) + labs(x = "Samples", y = "Cells",fontface = "bold") + 
scale_color_manual(values = c("#c90000",  "#ea9800",  "#04a327","#0073C2FF", "#6e0087")) + scale_fill_manual(values = c("#c90000",  "#ea9800",  "#04a327","#0073C2FF", "#6e0087"))
dev.off()

df <- df[df$Cell.Type!="Putative Cancer cells",]
df_full <- df[df$variable!="PE23",]

pdf("DGIST_full.pdf", width=10, height=10)
p <- ggplot(df_full, aes(fill=Cell.Type, y=value, x=variable))+geom_bar(stat="identity")+geom_text(aes(label = percentage),size =4.5, position = position_stack(vjust = 0.5),colour="white",fontface = "bold")
base_size <- 20
p + ggtitle("DGIST # of Cells by Cell Types") + 
theme_grey(base_size = base_size) + labs(x = "Samples", y = "Cells",fontface = "bold") + 
scale_color_manual(values = c("#c90000",  "#ea9800",  "#04a327", "#6e0087")) + scale_fill_manual(values = c("#c90000",  "#ea9800",  "#04a327", "#6e0087"))
dev.off()



library(ggcorrplot)
df <- read.table("cor.txt", sep="\t", header=TRUE)
rownames(df) <- df$X
df <- df[-1]

corr <- round(cor(df),2)
cor.test(df$DGIST, df$CIBERSORT)

# Pearson's product-moment correlation

#data:  df$DGIST and df$CIBERSORT
#t = 7.6758, df = 13, p-value = 3.507e-06
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7324729 0.9683825
#sample estimates:
#      cor 
#0.9051184 

pdf("corrplot.pdf", width=10, height=10)
corrplot(corr, method = "pie",type="upper",tl.col="black", tl.srt=45,addCoef.col = "black")
dev.off()


df <- read.table("cor_Monocyte_cluster6.txt", sep="\t", header=TRUE)
rownames(df) <- df$X
df <- df[-1]
corr <- cor(df)

pdf("cor_Monocyte_cluster6.pdf", width=5, height=5)
title="Monocyte_cluster6"
corrplot(corr, method = "pie",type="upper",title=title,tl.col="black", tl.srt=45,addCoef.col = "black",tl.cex =0.7,mar=c(0,0,3,0))
dev.off()

#######################################################################

load("ensemblGenes2018-01-14.RData")
DGIST=readRDS(file="BC_seurat.rds")

#Tcell

TGgene='ENSG00000184271'

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

################################################
setwd("/storage2/Project/CSC/10X/DGIST_data")
load("/storage2/Project/CSC/10X/05_QC/PE2324252629/ensemblGenes2018-01-14.RData")
library(RColorBrewer)

library(Seurat)
library(scater)
library(dplyr)

DGIST=readRDS(file="BC_seurat.rds")


exp =  DGIST@scale.data["ENSG00000005471",]

#TGgene=CIBERSORT celltype marker gene 
TGgene=c('ENSG00000005471','ENSG00000150967','ENSG00000072818','ENSG00000087085','ENSG00000102575','ENSG00000042980','ENSG00000134028','ENSG00000156140','ENSG00000169252','ENSG00000204472','ENSG00000163568','ENSG00000161905','ENSG00000012779','ENSG00000116748','ENSG00000101280','ENSG00000164512','ENSG00000128383','ENSG00000239713','ENSG00000128284','ENSG00000221963','ENSG00000103569','ENSG00000128805','ENSG00000137486','ENSG00000141505','ENSG00000161944','ENSG00000104043','ENSG00000230223','ENSG00000172232','ENSG00000112182','ENSG00000153064','ENSG00000043039','ENSG00000127152','ENSG00000140379','ENSG00000110987','ENSG00000162373','ENSG00000125864','ENSG00000123095','ENSG00000023445','ENSG00000136573','ENSG00000138756','ENSG00000101425','ENSG00000157764','ENSG00000174672','ENSG00000109743','ENSG00000113303','ENSG00000173715','ENSG00000118292','ENSG00000171860','ENSG00000197405','ENSG00000134830','ENSG00000178538','ENSG00000164047','ENSG00000137757','ENSG00000150636','ENSG00000108702','ENSG00000181374','ENSG00000276409','ENSG00000102970','ENSG00000275385','ENSG00000172724','ENSG00000115009','ENSG00000102962','ENSG00000274736','ENSG00000275302','ENSG00000271503','ENSG00000108688','ENSG00000108700','ENSG00000118971','ENSG00000184451','ENSG00000121807','ENSG00000183625','ENSG00000160791','ENSG00000112486','ENSG00000126353','ENSG00000117281','ENSG00000134061','ENSG00000177455','ENSG00000158477','ENSG00000158485','ENSG00000158481','ENSG00000158473','ENSG00000158488','ENSG00000116824','ENSG00000090659','ENSG00000012124','ENSG00000122223','ENSG00000198821','ENSG00000139193','ENSG00000178562','ENSG00000167851','ENSG00000105383','ENSG00000104894','ENSG00000004468','ENSG00000167286','ENSG00000198851','ENSG00000160654','ENSG00000010610','ENSG00000101017',
  'ENSG00000102245','ENSG00000110448','ENSG00000013725','ENSG00000129226','ENSG00000110848','ENSG00000173762','ENSG00000125726','ENSG00000137101','ENSG00000105369','ENSG00000007312','ENSG00000121594','ENSG00000114013','ENSG00000153563','ENSG00000172116','ENSG00000153283','ENSG00000158825','ENSG00000164045','ENSG00000154162','ENSG00000148600','ENSG00000105810','ENSG00000170956','ENSG00000124469','ENSG00000205923','ENSG00000126759','ENSG00000133048','ENSG00000064886','ENSG00000182022','ENSG00000147119','ENSG00000105205','ENSG00000153923','ENSG00000132514','ENSG00000069493','ENSG00000111729','ENSG00000172243','ENSG00000155962','ENSG00000092009','ENSG00000171812','ENSG00000206561','ENSG00000163751','ENSG00000117322','ENSG00000146592','ENSG00000096006','ENSG00000109943','ENSG00000100122','ENSG00000184371','ENSG00000164400','ENSG00000119535','ENSG00000077984','ENSG00000163599','ENSG00000100448','ENSG00000172543','ENSG00000169245','ENSG00000169248','ENSG00000156234','ENSG00000163734','ENSG00000163735','ENSG00000138755','ENSG00000163464','ENSG00000180871','ENSG00000160683','ENSG00000172215','ENSG00000147231','ENSG00000135929','ENSG00000111012','ENSG00000276644','ENSG00000035664','ENSG00000164935','ENSG00000164821','ENSG00000170456','ENSG00000100150','ENSG00000065357','ENSG00000278535','ENSG00000108771','ENSG00000167261','ENSG00000197635','ENSG00000134765','ENSG00000158050','ENSG00000145088','ENSG00000105246','ENSG00000184349','ENSG00000122877','ENSG00000197561','ENSG00000159023','ENSG00000072134','ENSG00000134954','ENSG00000117036','ENSG00000124019','ENSG00000185442','ENSG00000164125','ENSG00000117560','ENSG00000135722','ENSG00000179639','ENSG00000104921','ENSG00000072694','ENSG00000162747','ENSG00000085265','ENSG00000132704','ENSG00000182511',
  'ENSG00000126262','ENSG00000090554','ENSG00000119686','ENSG00000125740','ENSG00000049768','ENSG00000171051','ENSG00000171049','ENSG00000187474','ENSG00000111816','ENSG00000151474','ENSG00000126391','ENSG00000180340','ENSG00000104290','ENSG00000197093','ENSG00000166573','ENSG00000162676','ENSG00000099998','ENSG00000010310','ENSG00000176533','ENSG00000115523','ENSG00000076716','ENSG00000183671','ENSG00000174946','ENSG00000125245','ENSG00000169508','ENSG00000183150','ENSG00000170128','ENSG00000140030','ENSG00000100351','ENSG00000228315',
  'ENSG00000197465','ENSG00000145649','ENSG00000100453','ENSG00000100450','ENSG00000113088','ENSG00000197540','ENSG00000084110','ENSG00000101336','ENSG00000140287','ENSG00000163666','ENSG00000152804','ENSG00000177374','ENSG00000277075','ENSG00000273802','ENSG00000160883','ENSG00000241106','ENSG00000196735','ENSG00000213652','ENSG00000150540','ENSG00000105991','ENSG00000163106','ENSG00000173083','ENSG00000196639','ENSG00000173110','ENSG00000135914','ENSG00000003147','ENSG00000163600','ENSG00000131203','ENSG00000137959','ENSG00000186803','ENSG00000111537','ENSG00000211898','ENSG00000211891','ENSG00000211899','ENSG00000211592','ENSG00000206066','ENSG00000140749','ENSG00000113302','ENSG00000081985','ENSG00000112115','ENSG00000115604','ENSG00000115607','ENSG00000115008','ENSG00000125538','ENSG00000115602','ENSG00000138684','ENSG00000111536','ENSG00000134460','ENSG00000100385','ENSG00000164399','ENSG00000113520','ENSG00000077238','ENSG00000113525','ENSG00000091181','ENSG00000104432','ENSG00000168685','ENSG00000145839','ENSG00000140968','ENSG00000113263','ENSG00000177272','ENSG00000178342','ENSG00000125498','ENSG00000189013','ENSG00000221957','ENSG00000240403','ENSG00000111796','ENSG00000205810',
  'ENSG00000183542','ENSG00000134539','ENSG00000150045','ENSG00000139187','ENSG00000213809','ENSG00000219941','ENSG00000115919','ENSG00000089692','ENSG00000167618','ENSG00000078081','ENSG00000213658','ENSG00000182866','ENSG00000138795','ENSG00000138039','ENSG00000239998',
  'ENSG00000239961','ENSG00000131042','ENSG00000203896','ENSG00000281005','ENSG00000118308','ENSG00000204482','ENSG00000226979','ENSG00000227507','ENSG00000213316','ENSG00000112799','ENSG00000122224','ENSG00000185247','ENSG00000111837','ENSG00000111885','ENSG00000172469','ENSG00000073803','ENSG00000104814','ENSG00000168067','ENSG00000164114','ENSG00000173926','ENSG00000019169','ENSG00000105613','ENSG00000165471','ENSG00000103313','ENSG00000112818','ENSG00000257335','ENSG00000243156','ENSG00000262406','ENSG00000008516','ENSG00000100985','ENSG00000163563','ENSG00000184313','ENSG00000156738','ENSG00000149534','ENSG00000149516','ENSG00000110077','ENSG00000178860','ENSG00000059728','ENSG00000118513','ENSG00000170476','ENSG00000168060','ENSG00000116701','ENSG00000204475','ENSG00000123405','ENSG00000165028','ENSG00000105374','ENSG00000162711','ENSG00000135577','ENSG00000086288','ENSG00000167207','ENSG00000074771','ENSG00000130751','ENSG00000196436','ENSG00000135838','ENSG00000119508','ENSG00000162068','ENSG00000198400','ENSG00000085840','ENSG00000099985','ENSG00000108405','ENSG00000083454','ENSG00000078589','ENSG00000181631','ENSG00000174944','ENSG00000175591','ENSG00000159339','ENSG00000137819','ENSG00000115687','ENSG00000009709','ENSG00000163346','ENSG00000204965','ENSG00000188389','ENSG00000197646','ENSG00000095464','ENSG00000152256','ENSG00000008438','ENSG00000100100','ENSG00000078795','ENSG00000144837','ENSG00000146070','ENSG00000149527','ENSG00000166289',
  'ENSG00000126822','ENSG00000183395','ENSG00000168081','ENSG00000163736','ENSG00000110841','ENSG00000180644','ENSG00000186652','ENSG00000135362','ENSG00000242221','ENSG00000168229','ENSG00000125384','ENSG00000160013','ENSG00000213402','ENSG00000144724','ENSG00000213413','ENSG00000115828','ENSG00000041353','ENSG00000116191','ENSG00000185989','ENSG00000068831','ENSG00000152689','ENSG00000107551','ENSG00000117602','ENSG00000143839','ENSG00000102032','ENSG00000169891','ENSG00000090104','ENSG00000127074','ENSG00000169385','ENSG00000169413','ENSG00000165496','ENSG00000225093','ENSG00000052749','ENSG00000114767','ENSG00000134321','ENSG00000196218','ENSG00000163221','ENSG00000180739','ENSG00000155307','ENSG00000169432','ENSG00000075826','ENSG00000188404','ENSG00000184702','ENSG00000164402','ENSG00000129158','ENSG00000183918','ENSG00000088827','ENSG00000142178','ENSG00000089012','ENSG00000137078','ENSG00000154839','ENSG00000141293','ENSG00000117090','ENSG00000158714','ENSG00000074803','ENSG00000221955','ENSG00000110446','ENSG00000160326','ENSG00000130876','ENSG00000137571','ENSG00000103056','ENSG00000130768','ENSG00000185338','ENSG00000079263','ENSG00000061656','ENSG00000269404','ENSG00000107742','ENSG00000126752','ENSG00000064225','ENSG00000136840','ENSG00000111728','ENSG00000035720','ENSG00000127954','ENSG00000168952','ENSG00000233402','ENSG00000073861','ENSG00000081059','ENSG00000100721','ENSG00000135605','ENSG00000129566','ENSG00000104055','ENSG00000137462','ENSG00000196664','ENSG00000101916','ENSG00000121895','ENSG00000125355','ENSG00000123610','ENSG00000173535','ENSG00000141655','ENSG00000240505','ENSG00000048462','ENSG00000186827','ENSG00000125735','ENSG00000050730','ENSG00000172236','ENSG00000277734','ENSG00000076604','ENSG00000163519',
  'ENSG00000211789','ENSG00000211788','ENSG00000211791','ENSG00000211801','ENSG00000211795','ENSG00000211793','ENSG00000211751','ENSG00000211829','ENSG00000124731','ENSG00000095970','ENSG00000112195','ENSG00000071575','ENSG00000130529','ENSG00000119121','ENSG00000165409','ENSG00000075234','ENSG00000074966','ENSG00000077498','ENSG00000160185','ENSG00000242366','ENSG00000197888','ENSG00000100373','ENSG00000136059','ENSG00000112299','ENSG00000112303','ENSG00000093134','ENSG00000128218','ENSG00000111186','ENSG00000154764','ENSG00000115085','ENSG00000124256','ENSG00000205189','ENSG00000011590','ENSG00000152518','ENSG00000176293','ENSG00000197279','ENSG00000204789','ENSG00000159885','ENSG00000187607','ENSG00000083812','ENSG00000198342')

for (i in TGgene){
  test = DGIST@scale.data[i,]
  exp <- rbind(exp,test)
}
exp_real <-exp[rownames(exp)!="exp",]
rownames(exp_real) <- TGgene

dim(exp_real)
[1]   529 55788
write.table(exp_real, file = "DGIST_CIBERSORT.txt", sep="\t", row.names=TRUE)
exp_real <- read.table(file = "DGIST_CIBERSORT.txt", sep="\t", header=TRUE)

CIBERSORT <- read.table(file="CIBERSORT_marker_gene_compare_DGIST.txt", sep='\t', header=TRUE)
rownames(CIBERSORT) <- CIBERSORT[,1]
CIBERSORT <- CIBERSORT[,-1]
> CIBERSORT$
CIBERSORT$Gene.symbol                   CIBERSORT$NK.cells.activated
CIBERSORT$B.cells.naive                 CIBERSORT$Monocytes
CIBERSORT$B.cells.memory                CIBERSORT$Macrophages.M0
CIBERSORT$Plasma.cells                  CIBERSORT$Macrophages.M1
CIBERSORT$T.cells.CD8                   CIBERSORT$Macrophages.M2
CIBERSORT$T.cells.CD4.naive             CIBERSORT$Dendritic.cells.resting
CIBERSORT$T.cells.CD4.memory.resting    CIBERSORT$Dendritic.cells.activated
CIBERSORT$T.cells.CD4.memory.activated  CIBERSORT$Mast.cells.resting
CIBERSORT$T.cells.follicular.helper     CIBERSORT$Mast.cells.activated
CIBERSORT$T.cells.regulatory..Tregs.    CIBERSORT$Eosinophils
CIBERSORT$T.cells.gamma.delta           CIBERSORT$Neutrophils
CIBERSORT$NK.cells.resting              

#naiveB <- as.matrix(CIBERSORT$B.cells.naive)
#rownames(naiveB) <- rownames(CIBERSORT)

Neut <- as.matrix(CIBERSORT$Neutrophils)
rownames(Neut) <- rownames(CIBERSORT)


#naiveB_corr = c()
#memoryB_corr = c()
#plasma_corr = c()
#restingMemoryCD4T_corr = c()
#activatedMemoryCD4T_corr = c()
#folliTh_corr = c()
#Treg_corr = c()
#restingNK_corr = c()
Neut_corr = c()

for (i in c(1:length(colnames(exp_real)))){
table <- cbind(exp_real[,i],Neut[,1])
colnames(table) <- c("DGIST","CIBERSORT")
corr <- cor(table)[1,2]
Neut_corr <- cbind(Neut_corr,corr)
}
colnames(Neut_corr) <- colnames(exp_real)

celltypes = "Neutrophils"

pdf(paste0("DGIST_",celltypes,".pdf"),width=12, height=9)
df = data.frame(x=DGIST@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=DGIST@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=Neut_corr[1,])
ggplot(df,aes(x=x, y=y, colour=expression),mar=c(0,0,3,0)) + 
  ggtitle(paste0(celltypes, " Correlation")) +
  geom_point(size=0.5) + 
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



#DGIST vs Xcell : gene exp / cell type 
TGgene=c('ALDH1A2','ALOX15','CCL13','CCL17','CCL18','CCL23','CCL24','CD209','CD86','CLEC10A','F13A1','FCER2','IL3RA','SPINT2')

for (i in TGgene){
  print(ensemblGenes[ensemblGenes$external_gene_name == i,'ensembl_gene_id'])
}


exp =  DGIST@scale.data["ENSG00000128918",]

#TGgene=CIBERSORT celltype marker gene 
TGgene=c("ENSG00000161905","ENSG00000181374","ENSG00000102970","ENSG00000275385","ENSG00000274736","ENSG00000106178","ENSG00000090659","ENSG00000114013","ENSG00000132514","ENSG00000124491","ENSG00000104921","ENSG00000185291","ENSG00000167642")

for (i in TGgene){
  test = DGIST@scale.data[i,]
  exp <- rbind(exp,test)
}

celltypes = "iDC"

pdf(paste0("DGIST_Xcell_",celltypes,".pdf"),width=12, height=9)
df = data.frame(x=DGIST@dr$tsne@cell.embeddings[, "tSNE_1"], 
                y=DGIST@dr$tsne@cell.embeddings[, "tSNE_2"], 
                expression=colSums(exp))
ggplot(df,aes(x=x, y=y, colour=expression),mar=c(0,0,3,0)) + 
  ggtitle(paste0(celltypes, " gene exp")) +
  geom_point(size=0.5) + 
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

#####################################################
#####################################################
immune cell pop 따로
#####################################################
setwd('/storage2/Project/CSC/10X/DGIST_data')
DGIST=readRDS(file="BC_seurat.rds")

data_to_write_out <- as.data.frame(as.matrix(DGIST@scale.data))

fwrite(x = data_to_write_out, row.names = TRUE, file = "outfile.csv")

write.csv(BC_seurat.markers.genesymbol, file="PE26_BC_5000_seurat.markers.genesymbol.csv")

