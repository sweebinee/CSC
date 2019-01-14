# samplelist <- read.delim("C:/Users/jaewon/Desktop/2017work/PDX/sample_trio72list.txt", header = T, stringsAsFactors = F, sep="\t", na.strings = "")
# tumorlist <- samplelist$tumor[!is.na(samplelist$tumor)]
# xenolist <- samplelist[,c("PDX0","PDX1","PDX2","PDX3")][!is.na(samplelist[,c("PDX0","PDX1","PDX2","PDX3")])]
# allsamplelist <- c(tumorlist, xenolist)

samplelist = c('PE24','PE25','PE26','PE29','PE32','PE36')
setwd("/storage2/Project/CSC/WES/04_CNV/EXCAVATOR2/W_60000")

chrLen <- read.delim("/storage2/Project/CSC/WES/04_CNV/EXCAVATOR2/chrLen.txt", header = F, stringsAsFactors = F)
chrLen$V1 <- gsub("chr","",chrLen$V1) 
chrLen$V3 <- 0 
for(i in 2:nrow(chrLen)){
  chrLen[1,"V3"] <- chrLen[1,"V2"]
  chrLen[i,"V3"] <- sum(as.numeric(chrLen$V2[1:i]))
}
chrLen$idx <- 0
for(i in 2:nrow(chrLen)){
  chrLen[1,"idx"] <- chrLen[1,"V3"]/2
  chrLen[i,"idx"] <- (chrLen[i,"V3"]+chrLen[i-1,"V3"])/2
} 


for(i in 1:nrow(samplelist)){

  tsample <- samplelist[i,"tumor"]
  xsample <- samplelist[i,"PDX0"]
  
  setwd(paste0("C:/Users/jaewon/Desktop/2017work/PDX/WES/Excavator2/Result/",tsample))
  tt <- read.table(paste0("HSLMResults_",tsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  zz <- read.table(paste0("FastCallResults_",tsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  
  setwd(paste0("C:/Users/jaewon/Desktop/2017work/PDX/WES/Excavator2/Result/",xsample))
  xx <- read.table(paste0("HSLMResults_",xsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  ww <- read.table(paste0("FastCallResults_",xsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  
  
  for(chr in 2:(nrow(chrLen)-1)){
    tt[which(tt$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")] <- tt[which(tt$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")]+chrLen[chr-1,"V3"]
    xx[which(xx$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")] <- xx[which(xx$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")]+chrLen[chr-1,"V3"]
  }
  
  log2RSeqT <- as.numeric(as.character(tt[,5]))
  chrSeqT <- as.character(tt[,1])
  SegSeqT <- as.numeric(as.character(tt[,6]))
  PositionSeqT <- as.numeric(as.character(tt[,2]))
  TargetSeqT <- as.character(as.character(tt[,7]))
  
  log2RSeqX <- as.numeric(as.character(xx[,5]))
  chrSeqX <- as.character(xx[,1])
  SegSeqX <- as.numeric(as.character(xx[,6]))
  PositionSeqX <- as.numeric(as.character(xx[,2]))
  TargetSeqX <- as.character(as.character(xx[,7]))
  
  
  for(chr in 2:(nrow(chrLen)-1)){
    zz[which(zz$Chromosome==chrLen[chr,"V1"]),c("Start","End")] <- zz[which(zz$Chromosome==chrLen[chr,"V1"]),c("Start","End")]+chrLen[chr-1,"V3"]
    ww[which(ww$Chromosome==chrLen[chr,"V1"]),c("Start","End")] <- ww[which(ww$Chromosome==chrLen[chr,"V1"]),c("Start","End")]+chrLen[chr-1,"V3"]
  }
 
  chrCallT <- as.character(zz[,1])
  StartCallT <- as.numeric(as.character(zz[,2]))
  EndCallT <- as.numeric(as.character(zz[,3]))
  CallT <- as.numeric(as.character(zz[,7]))
  
  indINt <- which(TargetSeqT=="IN")
  indOUTt <- which(TargetSeqT=="OUT")

  
  chrCallX <- as.character(ww[,1])
  StartCallX <- as.numeric(as.character(ww[,2]))
  EndCallX <- as.numeric(as.character(ww[,3]))
  CallX <- as.numeric(as.character(ww[,7]))
  
  indINx <- which(TargetSeqX=="IN")
  indOUTx <- which(TargetSeqX=="OUT")
  
#  ------------------------------------------------------------------------
  
  
  png(filename=paste0("C:/Users/jaewon/Desktop/2017work/PDX/WES/Excavator2/Plot/",tsample,"_black.png"), width = 2000, height = 800, res=75)
  par(mfrow=c(2,1))
  
  if (length(indOUTt)!=0)
  {
    PositionSeqINt <- PositionSeqT[indINt]
    PositionSeqOUTt <- PositionSeqT[indOUTt]
    log2RSeqINt <- log2RSeqT[indINt]
    log2RSeqOUTt <- log2RSeqT[indOUTt]
    SegSeqINt <- SegSeqT[indINt]
    SegSeqOUTt <- SegSeqT[indOUTt]
    
    
    plot(PositionSeqOUTt,log2RSeqOUTt,ylim=c(-3,3),main=tsample,pch=19,cex=0.3,xlab="Position",ylab="log2ratio",col="blue",xaxt='n')
    points(PositionSeqINt,log2RSeqINt,col="black",pch=19,cex=0.3)
    lines(PositionSeqT, SegSeqT,lwd=2,col="red")
    abline(h=0,lty=2,lwd=1,col="green")
  }
  
  if (length(indOUTt)==0)
  {
    plot(PositionSeqT,log2RSeqT,ylim=c(-3,3),main=tsample,pch=19,cex=0.3,xlab="Position",ylab="log2ratio",col="blue",xaxt='n')
    lines(PositionSeqT, SegSeqT,lwd=2,col="red")
    abline(h=0,lty=2,lwd=1,col="green")
  }
  
  for(d in 1:23){
    abline(v=chrLen[d,"V3"]);
  }
  axis(1, at=chrLen$idx[1:23], labels=c(1:22,"X"),tick = F)
  
  
  ###
  
  if (length(indOUTx)!=0)
  {
    PositionSeqINx <- PositionSeqX[indINx]
    PositionSeqOUTx <- PositionSeqX[indOUTx]
    log2RSeqINx <- log2RSeqX[indINx]
    log2RSeqOUTx <- log2RSeqX[indOUTx]
    SegSeqINx <- SegSeqX[indINx]
    SegSeqOUTx <- SegSeqX[indOUTx]
    
    
    plot(PositionSeqOUTx,log2RSeqOUTx,ylim=c(-3,3),main=xsample,pch=19,cex=0.3,xlab="Position",ylab="log2ratio",col="blue",xaxt='n')
    points(PositionSeqINx,log2RSeqINx,col="black",pch=19,cex=0.3)
    lines(PositionSeqX, SegSeqX,lwd=2,col="red")
    abline(h=0,lty=2,lwd=1,col="green")
  }
  
  if (length(indOUTx)==0)
  {
    plot(PositionSeqX,log2RSeqX,ylim=c(-3,3),main=xsample,pch=19,cex=0.3,xlab="Position",ylab="log2ratio",col="blue",xaxt='n')
    lines(PositionSeqX, SegSeqX,lwd=2,col="red")
    abline(h=0,lty=2,lwd=1,col="green")
  }
  
  for(d in 1:23){
    abline(v=chrLen[d,"V3"]);
  }
  axis(1, at=chrLen$idx[1:23], labels=c(1:22,"X"),tick = F)
  
  
  dev.off()
  
}
  
  

#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------


for(i in 1:length(samplelist)){
  tsample <- samplelist[i]
  #xsample <- samplelist[i,"PDX0"]
  
  setwd(paste0("/storage2/Project/CSC/WES/04_CNV/EXCAVATOR2/",tsample,"_EXCAVATOR2/Results/",tsample))
  tt <- read.table(paste0("HSLMResults_",tsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  zz <- read.table(paste0("FastCallResults_",tsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  
  #setwd(paste0("C:/Users/jaewon/Desktop/2017work/PDX/WES/Excavator2/Result/",xsample))
  #xx <- read.table(paste0("HSLMResults_",xsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  #ww <- read.table(paste0("FastCallResults_",xsample,".txt"),sep="\t",quote="\"",fill=T,header=T, stringsAsFactors = F)
  
  for(chr in 2:(nrow(chrLen)-1)){
    tt[which(tt$Chromosome==paste0("chr",chrLen[chr,"V1"])),c("Position","Start","End")] <- tt[which(tt$Chromosome==paste0("chr",chrLen[chr,"V1"])),c("Position","Start","End")]+chrLen[chr-1,"V3"]
  #  xx[which(xx$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")] <- xx[which(xx$Chromosome==chrLen[chr,"V1"]),c("Position","Start","End")]+chrLen[chr-1,"V3"]
  }
  
  log2RSeqT <- as.numeric(as.character(tt[,5]))
  chrSeqT <- as.character(tt[,1])
  SegSeqT <- as.numeric(as.character(tt[,6]))
  PositionSeqT <- as.numeric(as.character(tt[,2]))
  TargetSeqT <- as.character(as.character(tt[,7]))
  
  #log2RSeqX <- as.numeric(as.character(xx[,5]))
  #chrSeqX <- as.character(xx[,1])
  #SegSeqX <- as.numeric(as.character(xx[,6]))
  #PositionSeqX <- as.numeric(as.character(xx[,2]))
  #TargetSeqX <- as.character(as.character(xx[,7]))
  
  for(chr in 2:(nrow(chrLen)-1)){
    zz[which(zz$Chromosome==paste0("chr",chrLen[chr,"V1"])),c("Start","End")] <- zz[which(zz$Chromosome==paste0("chr",chrLen[chr,"V1"])),c("Start","End")]+chrLen[chr-1,"V3"]
  #  ww[which(ww$Chromosome==chrLen[chr,"V1"]),c("Start","End")] <- ww[which(ww$Chromosome==chrLen[chr,"V1"]),c("Start","End")]+chrLen[chr-1,"V3"]
  }
  
  indINt <- which(TargetSeqT=="IN")
  indOUTt <- which(TargetSeqT=="OUT")
  
  #indINx <- which(TargetSeqX=="IN")
  #indOUTx <- which(TargetSeqX=="OUT")
  
  #  ------------------------------------------------------------------------
  
  
  png(filename=paste0("/storage2/Project/CSC/WES/04_CNV/EXCAVATOR2/W_10K_",tsample,"_grey.png"), width = 2000, height = 800, res=75)
#  par(mfrow=c(2,1))
  
  plot(PositionSeqT,log2RSeqT,ylim=c(-2,2),main=tsample,pch=19,cex=0.3,lwd=0.1,col="white",xlab="chromosome",ylab="log2ratio",xaxt='n')

  chrCallT <- as.character(zz[,1])
  StartCallT <- as.numeric(as.character(zz[,2]))
  EndCallT <- as.numeric(as.character(zz[,3]))
  CallT <- as.numeric(as.character(zz[,7]))
  
  for (j in 1:nrow(zz)){
    if (CallT[j]==1){
      rect(xleft=StartCallT[j], ybottom=0, xright=EndCallT[j], ytop=1, density = NA, angle = 45, col = "red",border="red")
    }
    if (CallT[j]==2){
      rect(xleft=StartCallT[j], ybottom=0, xright=EndCallT[j], ytop=2, density = NA, angle = 45, col = "darkred",border="darkred")
    }
    if (CallT[j]==-1){
      rect(xleft=StartCallT[j], ybottom=-1, xright=EndCallT[j], ytop=0, density = NA, angle = 45, col = "green",border="green")
    }
    if (CallT[j]==-2){
      rect(xleft=StartCallT[j], ybottom=-2, xright=EndCallT[j], ytop=0, density = NA, angle = 45, col = "darkgreen",border="darkgreen")
    }
  }
  
  #points(PositionSeqT,log2RSeqT,ylim=c(-3,3),pch=19,cex=0.1,col="grey")
  abline(h=0,lwd=1,col="black")
  for(d in 1:23){
    abline(v=chrLen[d,"V3"]);
  }
  axis(1, at=chrLen$idx[1:23], labels=c(1:22,"X"),tick = F)
  
  dev.off()
  
}


#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------






