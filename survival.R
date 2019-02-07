cd /storage2/Project/CSC/10X/DGIST_data/01_survival
cut -f2,14,16,25 clinical.tsv > BRCA_clinical.txt
#submitter_id	vital_status	days_to_death	days_to_last_follow_up


setwd("/storage2/Project/CSC/10X/DGIST_data/01_survival")
library(survival)
library(survminer)

survival <- read.table("BRCA_clinical.txt",header=TRUE,stringsAsFactors=FALSE) #1098 patients 945 alive, 153 dead
survival$vital_status <- ifelse(survival$vital_status == "alive", '0', '1') # alive=0, dead=1 

for(i in 1:nrow(survival)){
	if(survival$vital_status[i]=='0'){
		survival$days_to_death[i] <- survival$days_to_last_follow_up[i]
	}
}
#survival[survival$days_to_death=="--",]
survival <- survival[!(survival$days_to_death=="--"),] #1095 patients 944 alive, 151 dead
#write.table(survival, file="BRCA_clinical_final.txt", sep="\t", col.names=TRUE, row.names=FALSE)
survival <- read.table("BRCA_clinical_final.txt",header=TRUE,stringsAsFactors=FALSE)

##download all patients' transcriptome profiling > gene exp quant > RNA-seq > FPKM data from gdc portal Repository
#gz.sh

src_dir='/storage2/Project/CSC/10X/DGIST_data/01_survival/TCGA_BRCA_trans_profiling/tmp'
src_file=list.files(src_dir)
name <- strsplit(src_file[1],split=".txt")[[1]]
exp<-read.table(paste(src_dir,"/",src_file[1],sep=""),col.names=c("gene",name),sep="\t",header=FALSE,stringsAsFactors=FALSE)
for(i in 2:length(src_file)){
	name <- strsplit(src_file[i],split=".txt")[[1]]
	data <- read.table(paste(src_dir,"/",src_file[i],sep=""),col.names=c("gene",name),sep="\t",header=FALSE,stringsAsFactors=FALSE)
	exp <- merge(x=exp, y=data, by='gene', all=TRUE)
	print(i)
}
rownames(exp) <- exp[,1]
exp <- exp[,-1]
for(i in 1:nrow(exp)){
	rownames(exp)[i] <- strsplit(rownames(exp)[i],split="\\.")[[1]][1]
}

#some txt files have sentences at the front|back of the files
#removed it by sed '/id/d', manually

cut -f1,6 gdc_sample_sheet.2019-01-14.tsv > sample.txt

#caseID 
sample <- read.table("/storage2/Project/CSC/10X/DGIST_data/01_survival/sample.txt",header=TRUE,stringsAsFactors=FALSE)
for(i in 1:length(colnames(exp))){
	file <- gsub("X","",gsub("\\.","-",colnames(exp)[i]))
	case <- sample[sample$FileID==file,"CaseID"]
	colnames(exp)[i] <- case
}
#write.table(exp,"BRCA_exp_ENSG_caseID.txt",sep="\t", col.names=TRUE, row.names=TRUE)
exp <- read.table("BRCA_exp_ENSG_caseID.txt",header=TRUE,sep="\t",row.names=1,stringsAsFactors=FALSE)

#remove duplicated caseID 
exp_rmdupCase<-exp[,-(grep("\\.[0-9]", colnames(exp)))] #1084 patients

load("/storage2/Project/CSC/10X/DGIST_data/ensemblGenes2018-11-16.RData")

target_list = c('MTATP6P1','MT-ND4L','MT-TP','MT-ND6','XIST','KCNQ1OT1','MTRNR2L1','MT-ND5','SPTBN1','AC027290.2','EIF3B','TRIM28','PRRC2A','RCC2','FLNB','MAP7D1','PTPRF','HUWE1','KHDC4','FAM120A','WEE1','FLNA','MT-ND3','FASN','AKAP13','RRBP1','MT-ND1','TFRC','DUS1L','LYPD3','BCYRN1','TSC22D2','AC023157.3','NEAT1','SLC7A5','CHERP','TAPBP','AFF4','IMPAD1','MT-RNR1','ALYREF','ATP2A2','MT-CYB','FAM84B','CAND1','COLGALT1','CDH1','MT-ND2','ARID1A','MT-ATP6','PABPN1','UBE2H','TENT5A','NKTR','KDM2A','NPTN','KIF1B','KMT2C','MT-CO1','MT-ND4','MT-CO2','NCKAP1','SLC38A2','MT-RNR2','SPTAN1','MAT2A','MT-CO3','DSG2','RHOB','CD81','IRF2BP2','MYH9','RYBP','TRPS1','CLIP1','RABL6','H1FX','PPP1R14B','MKNK2','DYNC1H1','ENAH','TCF7L1','EIF4G1','KMT2E','NFIC','CLTC','CTSZ','AKAP9','RAP2B','PITRM1','IRF1','CEBPD','ADAM15','BRD2','PTPRK','VEGFA','TMEM165','CCNL1','ANKRD54','ATP1B1','UBE2S','WSB1','PIM3','SF1','C6orf132','SON','PNISR','PPP1R15A','SUPT5H','IFRD1','GPRC5A','SIRT7','PRRC2C','COL6A1','PRKDC','ZFP36L1','SYNGR2','SLC52A2','SRM','HDGF','EDN1','VMP1','KPNB1','GNB1','AL118516.1','LSR','GADD45B','PRPF38B','C3','DSP','CYR61','HK2','TRA2B','CLDN4','MALAT1','CTSB','GIPC1','CAPNS1','RBM25','JUP','NOC2L','SOX4','SYNE2','SDC4','RBM39','SLC9A3R1','MLF2','EIF3A','XBP1','TMED9','EFNA1','AUP1','ACTN4','MTDH','KLF6','SDC1','ATP1A1','NDUFB9','HSPA5','CALM2','GAL','FUS','GNAS','HES1','WARS','GBP1','STMN1','GADD45A','JUN','IFI6','RPS6','YWHAZ','PEG10','HNRNPH1','ELF3','COX6A1','CALD1','RPL31','RPS3A','EEF1A1','MFGE8','RPLP0','GAPDH','TUBB','PTMA')
surv_p <- data.frame(gene=target_list)

for(target in target_list){
	ensg <- ensemblGenes[ensemblGenes$external_gene_name == target,'ensembl_gene_id']
	if(length(ensg) > 1){
		i = 1
	    while(!(ensg[i] %in% rownames(exp_rmdupCase))){
	      i <- i+ 1
	    }
	    ensg <- ensg[i]
	}
	target_exp <- as.data.frame(t(exp_rmdupCase[rownames(exp_rmdupCase)==ensg,]))
	q1 <- quantile(target_exp,probs=c(0.45,0.55),na.rm=TRUE)[1]
	q3 <- quantile(target_exp,probs=c(0.45,0.55),na.rm=TRUE)[2]

	high_group = gsub("\\.","-",rownames(subset(target_exp,target_exp>=q3)))
	low_group = gsub("\\.","-",rownames(subset(target_exp,target_exp<=q1)))
	surv_plot <- data.frame(patient=c(high_group,low_group),group=c(rep("high",length(high_group)),rep("low",length(low_group))))

	for(i in 1:nrow(surv_plot)){
		if(!(surv_plot$patient[i]%in%survival$submitter_id)){
			surv_plot[-i,]
			next
		}
		surv_plot$time[i] <- survival[survival$submitter_id==surv_plot$patient[i],"days_to_death"]
		surv_plot$status[i] <- survival[survival$submitter_id==surv_plot$patient[i],"vital_status"]
	}
	surv_plot$time <- as.numeric(surv_plot$time)
	surv_plot$status <- as.numeric(surv_plot$status)

	fit <- survfit(Surv(time=surv_plot$time,event=surv_plot$status)~group,data=surv_plot)
	surv_p[surv_p$gene==target,"pvalue"] <- surv_pvalue(fit)[2]
	#summary(fit)
	if(surv_pvalue(fit)[2]<=0.05){
		gg <- ggsurvplot(fit, data=surv_plot, palette=c("Red","Blue"), pval=TRUE, risk.table=T, title=paste0("BRCA : ",target," quantile 45%"))
		ggsave(file = paste0(target,"_Q45_surv_plt.png"), plot=print(gg))
	}
}


#########################################################################################


## expression boxplot
gene <- c(rep("CDH1",nrow(CDH1)),rep("SOX4",nrow(SOX4)),rep("KDM2A",nrow(KDM2A)))
expression <- c(CDH1$ENSG00000039068,SOX4$ENSG00000124766,KDM2A$ENSG00000173120)
exp_box <- data.frame(gene,expression) 

target="CDH1"
png(paste0(target,"_exp_box.png"))
ggplot(data=exp_box[gene==target,],aes(x=gene,y=expression,color=gene))+
geom_boxplot(outlier.alpha=0, color="black")+
geom_jitter(aes(x=gene),alpha=0.5)+theme_bw()+
labs(title="TCGA-BRCA ",target," expression distribution")
dev.off()

