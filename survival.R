cd /storage2/Project/CSC/10X/DGIST_data/01_survival
cut -f2,14,16,25 clinical.tsv > BRCA_clinical.txt
#submitter_id	vital_status	days_to_death	days_to_last_follow_up


setwd("/storage2/Project/CSC/10X/DGIST_data/01_survival")
library(survival)
library(survminer)

survival <- read.table("BRCA_clinical.txt",header=TRUE,stringsAsFactors=FALSE) #1098 patients 945 alive, 153 dead
survival$vital_status <- ifelse(survival$vital_status == "alive", '1', '0')

for(i in 1:nrow(survival)){
	if(survival$vital_status[i]=='1'){
		survival$days_to_death[i] <- survival$days_to_last_follow_up[i]
	}
}
#survival[survival$days_to_death=="--",]
survival <- survival[!(survival$days_to_death=="--"),] #1095 patients 944 alive, 151 dead
write.table(survival, file="BRCA_clinical_final.txt", sep="\t", col.names=TRUE, row.names=FALSE)

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



sample2$Expr = as.factor(sample2$Expr)
sample2$vital_status = as.numeric(sample2$vital_status)

  surv_object <- Surv(time = sample2$days_to_death, event = sample2$vital_status)
  fit1 <- survfit(surv_object ~ Expr, data = sample2)
  #summary(fit1)
  genename=colnames(sample2)[6]
  gg <- ggsurvplot(fit1, data = sample2, pval = TRUE,risk.table = T,title=paste0(" THCA ",genename))
  ggsave(file = paste0("./THCA/tumorexpr35/",genename,".jpg"), plot=print(gg))
