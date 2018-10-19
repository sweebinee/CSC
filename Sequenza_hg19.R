###installation
#R pakage
source("https://bioconductor.org/biocLite.R")
biocLite("sequenza")
#Sequenza-utils.py
git clone https://bitbucket.org/sequenza_tools/sequenza-utils
cd sequenza-utils
python setup.py test
python setup.py install

##Prepare inputs
#Alternatively, it is possible to use the output of VarScan2

###Generating pileup files form BAM files

###Generating a genome-wide GC content file
# is used to normalize the depth ratio

###Generating a seqz file
#A seqz file contains genotype information, alleles and mutation frequency, and other feature

#To reduce the size of the seqz file, we recommend the use of a binning function provided in sequenza-utils.py. This binning decreases the memory requirement to load the data into R, and it also speeds up the processing of the sample.
/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/01_mkseqz.sh


##Exploring the seqz file & depth ratio normalization details
#in R
setwd("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal")
library("sequenza")

PE17.data <- read.seqz("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE17.small.seqz.gz") #read all data at once
PE18.data <- read.seqz("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE18.small.seqz.gz") #read all data at once
PE20.data <- read.seqz("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE20.small.seqz.gz") #read all data at once
PE24.data <- read.seqz("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE24.small.seqz.gz") #read all data at once
#seqz.data <- read.seqz(PE17.small.seqz.gz, chr.name = "1") ## read chr1 only

str(PE17.data, vec.len = 2)
'data.frame':	4150956 obs. of  14 variables:
 $ chromosome     : chr  "chr1" "chr1" ...
 $ position       : int  10050 10055 10068 10074 10091 ...
 $ base.ref       : chr  "N" "T" ...
 $ depth.normal   : int  9 9 10 10 9 ...
 $ depth.tumor    : int  14 12 13 14 15 ...
 $ depth.ratio    : num  1.74 1.33 ...
 $ Af             : num  1 0.9 0.923 0.929 0.933 ...
 $ Bf             : num  0 0 0 0 0 ...
 $ zygosity.normal: chr  "hom" "hom" ...
 $ GC.percent     : num  50 50 50 50 50 ...
 $ good.reads     : num  47 10 13 14 15 ...
 $ AB.normal      : chr  "N" "T" ...
 $ AB.tumor       : chr  "." "C0.1" ...
 $ tumor.strand   : chr  "0" "C1.0" ...

###Quality Control
#Each aligned base, in the next generation sequencing, is associated with a quality score.
#The sequenza-utils.py software is capable of filtering out bases with a quality score lower then a specified value (default, 20). 
#The number of reads that have passed the filter is returned in the column good.reads, 
#while the depth.tumor column contains the raw depth indicated in the pileup.

###Normalization of depth ratio
#The GC content bias affects most of the samples;  however,  some samples are more biased than others.  
#We attempt to remove this bias by normalizing with the mean depth ratio value of a corresponding GC content value.

gc.stats <- gc.norm(x = PE17.data$depth.ratio, gc = PE17.data$GC.percent)
#gc.stats <- gc.sample.stats(PE17.data)
#str(gc.stats)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
PE17.data$adjusted.ratio <- PE17.data$depth.ratio / gc.vect[as.character(PE17.data$GC.percent)]

png("PE17_sequenza.png", width=1500, height=700 )
par(mfrow = c(1,2), cex = 1, las = 1, bty ='l')
matplot(gc.stats$gc.values, gc.stats$raw,type ='b', col = 1, pch = c(1, 19, 1), lty = c(2, 1, 2),xlab ='GC content (%)', ylab ='Uncorrected depth ratio')
legend('topright', legend = colnames(gc.stats$raw), pch = c(1, 19, 1))
hist2(PE17.data$depth.ratio, PE17.data$adjusted.ratio,breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),xlab ='Uncorrected depth ratio', ylab ='GC-adjusted depth ratio')
dev.off()

###Analyzing sequencing data with sequenza
#Extract the information from the .seqz file
$ zcat PE17.small.seqz.gz | sed '/random/d' | sed '/alt/d' | gzip > PE17.txt.gz

test <- sequenza.extract("PE17.txt.gz")
names(test)
[1] "BAF"         "ratio"       "mutations"   "segments"    "chromosomes" "gc"          "avg.depth"  

# Plot chromosome view with mutations, BAF, depth ratio and segments
for (i in c(1:23)){
png(paste0("PE17_",test$chromosomes[i],".png"), width=450, height=400 )
chromosome.view(mut.tab = test$mutations[[i]], baf.windows = test$BAF[[i]],
                 ratio.windows = test$ratio[[i]], min.N.ratio = 1,
                 segments = test$segments[[i]], main = test$chromosomes[i])
dev.off()}
#the calculated B allele frequency and depth  ratio  of  the  obtained  segments
CP.example <- sequenza.fit(test)
#Results of model fitting
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = "PE17", out.dir="PE17")

cint <- get.ci(CP.example)

png("PE17_CP.png")
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE, likThresh = c(0.95))
dev.off()

png("PE17_log_prob_ploidy.png")
par(mfrow = c(2,2))
cp.plot(CP.example)
cp.plot.contours(CP.example, add = TRUE)
plot(cint$values.cellularity, ylab = "Cellularity",xlab = "posterior probability", type = "n")
select <- cint$confint.cellularity[1] <= cint$values.cellularity[,2] &cint$values.cellularity[,2] <= cint$confint.cellularity[2]
polygon(y = c(cint$confint.cellularity[1], cint$values.cellularity[select, 2], cint$confint.cellularity[2]),x = c(0, cint$values.cellularity[select, 1], 0), col='red', border=NA)
lines(cint$values.cellularity)
abline(h = cint$max.cellularity, lty = 2, lwd = 0.5)
plot(cint$values.ploidy, xlab = "Ploidy",ylab = "posterior probability", type = "n")
select <- cint$confint.ploidy[1] <= cint$values.ploidy[,1] & cint$values.ploidy[,1] <= cint$confint.ploidy[2]
polygon(x = c(cint$confint.ploidy[1], cint$values.ploidy[select, 1], cint$confint.ploidy[2]),y = c(0, cint$values.ploidy[select, 2], 0), col='red', border=NA)
lines(cint$values.ploidy)
abline(v = cint$max.ploidy, lty = 2, lwd = 0.5)
dev.off()

# Call CNVs and mutations
cellularity <- cint$max.cellularity
cellularity
[1] 0.35

ploidy <- cint$max.ploidy
ploidy
[1] 2.4

avg.depth.ratio <- mean(test$gc$adj[, 2])
avg.depth.ratio
[1] 1

# Detect variant alleles (mutations)
mut.tab <- na.exclude(do.call(rbind, test$mutations))
mut.alleles <- mufreq.bayes(mufreq = mut.tab$F,depth.ratio = mut.tab$adjusted.ratio,cellularity = cellularity, ploidy = ploidy,avg.depth.ratio = avg.depth.ratio)
head(mut.alleles)
   CNn CNt Mt        LPP
7    2   3  1  -8.592913
4    2   2  1 -10.356419
41   2   2  1  -9.351969
42   2   2  1  -8.399611
71   2   3  1  -8.547686
72   2   3  1  -8.808447

head(cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")],mut.alleles))

# Detect copy number variations
seg.tab     <- na.exclude(do.call(rbind, test$segments))
cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)
head(cn.alleles)

seg.tab <- cbind(seg.tab, cn.alleles)
head(seg.tab)

for (i in c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")){
png(paste0("PE17_",test$chromosomes[i],"_alSpec.png"), width=450, height=400 )
chromosome.view(mut.tab = test$mutations[[i]], baf.windows = test$BAF[[i]],
                 ratio.windows = test$ratio[[i]],  min.N.ratio = 1,
                 segments = seg.tab[seg.tab$chromosome == i,],
                 main = test$chromosomes[i],
                 cellularity = cellularity, ploidy = ploidy,
                 avg.depth.ratio = avg.depth.ratio)
dev.off()}

seg_order <- seg.tab[order(chromosome),]
png("PE17_CNV_genome.png", width=700, height=250 )
genome.view(seg.cn = seg.tab, info.type = "CNt")
legend("bottomright", bty="n", c("Tumor copy number"),col = c("red"), inset = c(0, -0.4), pch=15, xpd = TRUE)
dev.off()

png("PE17_CNV_alSpec.png")
genome.view(seg.cn = seg.tab, info.type = "AB")
legend("bottomright", bty = "n", c("A-allele","B-allele"), col= c("red", "blue"),
        inset = c(0, -0.45), pch = 15, xpd = TRUE)
dev.off()

##################################################################################
#for pyclone input 
#############################
setwd("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal")
library("sequenza")

sequenza:::sequenza2PyClone

sequenza2PyClone <- function (mut.tab, seg.cn, sample.id, norm.cn = 2) {
    mut.tab <- cbind(mut.tab[, c("chromosome", "position", "good.reads", "F", "mutation")], CNt = NA, A = NA, B = NA)
    for (i in 1:nrow(seg.cn)) {
        pos.filt <- mut.tab$chromosome == seg.cn$chromosome[i] & 
            mut.tab$position >= seg.cn$start.pos[i] & mut.tab$position <= seg.cn$end.pos[i]
        mut.tab[pos.filt, c("CNt", "A", "B")] <- seg.cn[i, c("CNt", "A", "B")]
    }
    id <- paste(mut.tab$chromosome, mut.tab$position, sep = ":")
    var.counts <- round(mut.tab$good.reads * mut.tab$F, 0)
    nor.counts <- mut.tab$good.reads - var.counts
    pyclone.tsv <- data.frame(mutation_id = id, ref_counts = nor.counts, 
        var_counts = var.counts, normal_cn = norm.cn, minor_cn = mut.tab$B, 
        major_cn = mut.tab$A, variant_case = sample.id, variant_freq = mut.tab$F, 
        genotype = mut.tab$mutation)
    na.exclude(pyclone.tsv)
}


PE17.data <- read.seqz("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE17.small.seqz.gz") #read all data at once

gc.stats <- gc.norm(x = PE17.data$depth.ratio, gc = PE17.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
PE17.data$adjusted.ratio <- PE17.data$depth.ratio / gc.vect[as.character(PE17.data$GC.percent)]

$ zcat BC24.small.seqz.gz | sed '/_/d' | gzip > BC24.txt.gz

test <- sequenza.extract("/storage2/Project/CSC/WES/ver.hg19/03_longitudinal/PE17.txt.gz")

CP.example <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP.example, sample.id = "PE17", out.dir="PE17")
cint <- get.ci(CP.example)

cellularity <- cint$max.cellularity
ploidy <- cint$max.ploidy
avg.depth.ratio <- mean(test$gc$adj[, 2])

mut.tab <- na.exclude(do.call(rbind, test$mutations))
mut.alleles <- mufreq.bayes(mufreq = mut.tab$F,depth.ratio = mut.tab$adjusted.ratio,cellularity = cellularity, ploidy = ploidy,avg.depth.ratio = avg.depth.ratio)
head(cbind(mut.tab[,c("chromosome","position","F","adjusted.ratio", "mutation")],mut.alleles))

seg.tab <- na.exclude(do.call(rbind, test$segments))
cn.alleles <- baf.bayes(Bf = seg.tab$Bf, depth.ratio = seg.tab$depth.ratio, cellularity = cellularity, ploidy = ploidy, avg.depth.ratio = avg.depth.ratio)
head(cn.alleles)

seg.tab <- cbind(seg.tab, cn.alleles)
head(seg.tab)

result <- sequenza2PyClone(mut.tab, seg.cn=seg.tab, sample.id="PE17", norm.cn = 2)
write.table(result, file='PE17_broad_seq2pyC.tsv', quote=FALSE, sep='\t', col.names = NA)

#############################
#from TSV -> YAML mut
#############################
PyClone build_mutations_file PE17_broad_seq2pyC.tsv 