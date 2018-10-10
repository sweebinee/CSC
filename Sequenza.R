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
main_dir=/storage2/Project/CSC/WES/04_CNV/Sequenza
input_dir=/storage2/Project/CSC/WES/02_Preprocessing
ref_dir=/storage2/Project/source/ref_GRCh38

samtools mpileup −f $ref_dir/hg38.fa -Q 20 $input_dir/PE17_RECAL.bam | gzip > $main_dir/PE17.pileup.gz
samtools mpileup −f $ref_dir/hg38.fa -Q 20 $input_dir/PE18_RECAL.bam | gzip > $main_dir/PE18.pileup.gz
samtools mpileup −f $ref_dir/hg38.fa -Q 20 $input_dir/PE20_RECAL.bam | gzip > $main_dir/PE20.pileup.gz
samtools mpileup −f $ref_dir/hg38.fa -Q 20 $input_dir/PE24_RECAL.bam | gzip > $main_dir/PE24.pileup.gz
samtools mpileup −f $ref_dir/hg38.fa -Q 20 $input_dir/BC24_RECAL.bam | gzip > $main_dir/BC24.pileup.gz

###Generating a genome-wide GC content file
# is used to normalize the depth ratio
sequenza-utils gc_wiggle -w 50 -f $ref_dir/hg38.fa | gzip > $ref_dir/hg38.gc50Base.txt.gz

###Generating a seqz file
#A seqz file contains genotype information, alleles and mutation frequency, and other feature
sequenza-utils bam2seqz −gc $ref_dir/hg38.gc50Base.txt.gz \
-F $ref_dir/hg38.fa \
-p \
-n $main_dir/BC24.pileup.gz \
-t $main_dir/$sample".pileup.gz" \
-o $main_dir/$sample".seqz.gz" 
#To reduce the size of the seqz file, we recommend the use of a binning function provided in sequenza-utils.py. This binning decreases the memory requirement to load the data into R, and it also speeds up the processing of the sample.
sequenza-utils seqz_binning -w 50 -s $main_dir/$sample".seqz.gz" | gzip > $main_dir/$sample".small.seqz.gz"


##Exploring the seqz file & depth ratio normalization details
#in R
setwd("/storage2/Project/CSC/WES/04_CNV/Sequenza")
library("sequenza")

PE17.data <- read.seqz("/storage2/Project/CSC/WES/04_CNV/Sequenza/PE17.small.seqz.gz") #read all data at once
PE18.data <- read.seqz("/storage2/Project/CSC/WES/04_CNV/Sequenza/PE18.small.seqz.gz") #read all data at once
PE20.data <- read.seqz("/storage2/Project/CSC/WES/04_CNV/Sequenza/PE20.small.seqz.gz") #read all data at once
PE24.data <- read.seqz("/storage2/Project/CSC/WES/04_CNV/Sequenza/PE24.small.seqz.gz") #read all data at once
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
setwd("/storage2/Project/CSC/WES/04_CNV/Sequenza")
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

PE17.data <- read.seqz("/storage2/Project/CSC/WES/04_CNV/Sequenza/w_50/PE17.small.seqz.gz") #read all data at once

gc.stats <- gc.norm(x = PE17.data$depth.ratio, gc = PE17.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
PE17.data$adjusted.ratio <- PE17.data$depth.ratio / gc.vect[as.character(PE17.data$GC.percent)]

$ zcat PE17.small.seqz.gz | sed 's/_KI270721v1_random//g' | sed 's/_GL000009v2_random//g' | sed 's/_GL000194v1_random//g' | sed 's/_GL000225v1_random//g' | sed 's/_KI270722v1_random//g' | sed 's/_KI270723v1_random//g' | sed 's/_KI270724v1_random//g' | sed 's/_KI270725v1_random//g' | sed 's/_KI270726v1_random//g' | sed 's/_KI270727v1_random//g' | sed 's/_KI270728v1_random//g' | sed 's/_GL000205v2_random//g' | sed 's/_KI270729v1_random//g' | sed 's/_KI270706v1_random//g' | sed 's/_KI270708v1_random//g' | sed 's/_KI270709v1_random//g' | sed 's/_KI270711v1_random//g' | sed 's/_KI270712v1_random//g' | sed 's/_KI270713v1_random//g' | sed 's/_KI270714v1_random//g' | sed 's/_KI270731v1_random//g' | sed 's/_KI270732v1_random//g' | sed 's/_KI270733v1_random//g' | sed 's/_KI270734v1_random//g' | sed 's/_KI270735v1_random//g' | sed 's/_KI270736v1_random//g' | sed 's/_KI270737v1_random//g' | sed 's/_GL000221v1_random//g' | sed 's/_GL000008v2_random//g' | sed 's/_KI270717v1_random//g' | sed 's/_KI270718v1_random//g' | sed 's/_KI270719v1_random//g' | sed 's/_KI270720v1_random//g' |sed 's/_GL383545v1_alt//g' | sed 's/_GL383546v1_alt//g' | sed 's/_KI270824v1_alt//g' | sed 's/_KI270825v1_alt//g' | sed 's/_JH159136v1_alt//g' | sed 's/_JH159137v1_alt//g' | sed 's/_KI270829v1_alt//g' | sed 's/_KI270830v1_alt//g' | sed 's/_KI270831v1_alt//g' | sed 's/_KI270832v1_alt//g' | sed 's/_KI270902v1_alt//g' | sed 's/_KI270903v1_alt//g' | sed 's/_KI270927v1_alt//g' | sed 's/_GL383550v2_alt//g' | sed 's/_GL383553v2_alt//g' | sed 's/_GL877875v1_alt//g' | sed 's/_GL877876v1_alt//g' | sed 's/_KI270833v1_alt//g' | sed 's/_KI270834v1_alt//g' | sed 's/_KI270835v1_alt//g' | sed 's/_KI270837v1_alt//g' | sed 's/_KI270904v1_alt//g' | sed 's/_KI270838v1_alt//g' | sed 's/_KI270839v1_alt//g' | sed 's/_KI270840v1_alt//g' | sed 's/_KI270842v1_alt//g' | sed 's/_KI270844v1_alt//g' | sed 's/_KI270845v1_alt//g' | sed 's/_KI270846v1_alt//g' | sed 's/_KI270847v1_alt//g' | sed 's/_GL383554v1_alt//g' | sed 's/_GL383555v2_alt//g' | sed 's/_KI270848v1_alt//g' | sed 's/_KI270849v1_alt//g' | sed 's/_KI270850v1_alt//g' | sed 's/_KI270851v1_alt//g' | sed 's/_KI270852v1_alt//g' | sed 's/_KI270905v1_alt//g' | sed 's/_KI270906v1_alt//g' | sed 's/_GL383556v1_alt//g' | sed 's/_GL383557v1_alt//g' | sed 's/_KI270853v1_alt//g' | sed 's/_KI270854v1_alt//g' | sed 's/_KI270855v1_alt//g' | sed 's/_KI270856v1_alt//g' | sed 's/_GL000258v2_alt//g' | sed 's/_GL383563v3_alt//g' | sed 's/_GL383564v2_alt//g' | sed 's/_GL383566v1_alt//g' | sed 's/_JH159146v1_alt//g' | sed 's/_JH159147v1_alt//g' | sed 's/_JH159148v1_alt//g' | sed 's/_KI270857v1_alt//g' | sed 's/_KI270858v1_alt//g' | sed 's/_KI270859v1_alt//g' | sed 's/_KI270860v1_alt//g' | sed 's/_KI270861v1_alt//g' | sed 's/_KI270862v1_alt//g' | sed 's/_KI270907v1_alt//g' | sed 's/_KI270908v1_alt//g' | sed 's/_KI270909v1_alt//g' | sed 's/_KI270910v1_alt//g' | sed 's/_GL383567v1_alt//g' | sed 's/_GL383571v1_alt//g' | sed 's/_GL383572v1_alt//g' | sed 's/_KI270863v1_alt//g' | sed 's/_KI270911v1_alt//g' | sed 's/_GL000209v2_alt//g' | sed 's/_GL383573v1_alt//g' | sed 's/_GL383574v1_alt//g' | sed 's/_GL383575v2_alt//g' | sed 's/_GL383576v1_alt//g' | sed 's/_GL949746v1_alt//g' | sed 's/_GL949747v2_alt//g' | sed 's/_GL949748v2_alt//g' | sed 's/_GL949749v2_alt//g' | sed 's/_GL949750v2_alt//g' | sed 's/_GL949751v2_alt//g' | sed 's/_GL949752v1_alt//g' | sed 's/_GL949753v2_alt//g' | sed 's/_KI270865v1_alt//g' | sed 's/_KI270866v1_alt//g' | sed 's/_KI270867v1_alt//g' | sed 's/_KI270868v1_alt//g' | sed 's/_KI270882v1_alt//g' | sed 's/_KI270883v1_alt//g' | sed 's/_KI270884v1_alt//g' | sed 's/_KI270885v1_alt//g' | sed 's/_KI270886v1_alt//g' | sed 's/_KI270887v1_alt//g' | sed 's/_KI270890v1_alt//g' | sed 's/_KI270891v1_alt//g' | sed 's/_KI270915v1_alt//g' | sed 's/_KI270916v1_alt//g' | sed 's/_KI270917v1_alt//g' | sed 's/_KI270918v1_alt//g' | sed 's/_KI270920v1_alt//g' | sed 's/_KI270921v1_alt//g' | sed 's/_KI270922v1_alt//g' | sed 's/_KI270923v1_alt//g' | sed 's/_KI270929v1_alt//g' | sed 's/_KI270930v1_alt//g' | sed 's/_KI270938v1_alt//g' | sed 's/_GL383518v1_alt//g' | sed 's/_GL383519v1_alt//g' | sed 's/_GL383520v2_alt//g' | sed 's/_KI270759v1_alt//g' | sed 's/_KI270760v1_alt//g' | sed 's/_KI270761v1_alt//g' | sed 's/_KI270762v1_alt//g' | sed 's/_KI270763v1_alt//g' | sed 's/_KI270766v1_alt//g' | sed 's/_KI270892v1_alt//g' | sed 's/_KI270869v1_alt//g' | sed 's/_KI270870v1_alt//g' | sed 's/_KI270871v1_alt//g' | sed 's/_GL383579v2_alt//g' | sed 's/_GL383580v2_alt//g' | sed 's/_GL383581v2_alt//g' | sed 's/_KI270872v1_alt//g' | sed 's/_KI270873v1_alt//g' | sed 's/_GL383582v2_alt//g' | sed 's/_GL383583v2_alt//g' | sed 's/_KB663609v1_alt//g' | sed 's/_KI270875v1_alt//g' | sed 's/_KI270876v1_alt//g' | sed 's/_KI270877v1_alt//g' | sed 's/_KI270878v1_alt//g' | sed 's/_KI270879v1_alt//g' | sed 's/_KI270928v1_alt//g' | sed 's/_GL383521v1_alt//g' | sed 's/_GL383522v1_alt//g' | sed 's/_GL582966v2_alt//g' | sed 's/_KI270767v1_alt//g' | sed 's/_KI270768v1_alt//g' | sed 's/_KI270769v1_alt//g' | sed 's/_KI270772v1_alt//g' | sed 's/_KI270773v1_alt//g' | sed 's/_KI270774v1_alt//g' | sed 's/_KI270775v1_alt//g' | sed 's/_KI270776v1_alt//g' | sed 's/_KI270893v1_alt//g' | sed 's/_KI270894v1_alt//g' | sed 's/_GL383526v1_alt//g' | sed 's/_JH636055v2_alt//g' | sed 's/_KI270777v1_alt//g' | sed 's/_KI270779v1_alt//g' | sed 's/_KI270780v1_alt//g' | sed 's/_KI270781v1_alt//g' | sed 's/_KI270782v1_alt//g' | sed 's/_KI270784v1_alt//g' | sed 's/_KI270895v1_alt//g' | sed 's/_KI270924v1_alt//g' | sed 's/_KI270934v1_alt//g' | sed 's/_KI270935v1_alt//g' | sed 's/_KI270936v1_alt//g' | sed 's/_KI270937v1_alt//g' | sed 's/_GL000257v2_alt//g' | sed 's/_GL383527v1_alt//g' | sed 's/_KI270789v1_alt//g' | sed 's/_KI270790v1_alt//g' | sed 's/_KI270896v1_alt//g' | sed 's/_KI270925v1_alt//g' | sed 's/_GL339449v2_alt//g' | sed 's/_KI270791v1_alt//g' | sed 's/_KI270792v1_alt//g' | sed 's/_KI270793v1_alt//g' | sed 's/_KI270794v1_alt//g' | sed 's/_KI270795v1_alt//g' | sed 's/_KI270897v1_alt//g' | sed 's/_KI270898v1_alt//g' | sed 's/_GL000250v2_alt//g' | sed 's/_GL000251v2_alt//g' | sed 's/_GL000252v2_alt//g' | sed 's/_GL000253v2_alt//g' | sed 's/_GL000254v2_alt//g' | sed 's/_GL000255v2_alt//g' | sed 's/_GL000256v2_alt//g' | sed 's/_GL383533v1_alt//g' | sed 's/_KI270758v1_alt//g' | sed 's/_KI270797v1_alt//g' | sed 's/_KI270798v1_alt//g' | sed 's/_KI270800v1_alt//g' | sed 's/_KI270801v1_alt//g' | sed 's/_GL383534v2_alt//g' | sed 's/_KI270803v1_alt//g' | sed 's/_KI270804v1_alt//g' | sed 's/_KI270805v1_alt//g' | sed 's/_KI270806v1_alt//g' | sed 's/_KI270808v1_alt//g' | sed 's/_KI270809v1_alt//g' | sed 's/_KI270899v1_alt//g' | sed 's/_KI270811v1_alt//g' | sed 's/_KI270812v1_alt//g' | sed 's/_KI270813v1_alt//g' | sed 's/_KI270815v1_alt//g' | sed 's/_KI270816v1_alt//g' | sed 's/_KI270817v1_alt//g' | sed 's/_KI270818v1_alt//g' | sed 's/_KI270819v1_alt//g' | sed 's/_KI270820v1_alt//g' | sed 's/_KI270821v1_alt//g' | sed 's/_KI270822v1_alt//g' | sed 's/_KI270900v1_alt//g' | sed 's/_KI270926v1_alt//g' | sed 's/_GL383540v1_alt//g' | sed 's/_GL383541v1_alt//g' | sed 's/_GL383542v1_alt//g' | sed 's/_KI270823v1_alt//g' | sed 's/_KI270881v1_alt//g' | sed 's/_KI270913v1_alt//g' | gzip > PE17.txt.gz
$ zcat PE18.small.seqz.gz | sed 's/_KI270721v1_random//g' | sed 's/_GL000009v2_random//g' | sed 's/_GL000194v1_random//g' | sed 's/_GL000225v1_random//g' | sed 's/_KI270722v1_random//g' | sed 's/_KI270723v1_random//g' | sed 's/_KI270724v1_random//g' | sed 's/_KI270725v1_random//g' | sed 's/_KI270726v1_random//g' | sed 's/_KI270727v1_random//g' | sed 's/_KI270728v1_random//g' | sed 's/_GL000205v2_random//g' | sed 's/_KI270729v1_random//g' | sed 's/_KI270706v1_random//g' | sed 's/_KI270708v1_random//g' | sed 's/_KI270709v1_random//g' | sed 's/_KI270711v1_random//g' | sed 's/_KI270712v1_random//g' | sed 's/_KI270713v1_random//g' | sed 's/_KI270714v1_random//g' | sed 's/_KI270731v1_random//g' | sed 's/_KI270732v1_random//g' | sed 's/_KI270733v1_random//g' | sed 's/_KI270734v1_random//g' | sed 's/_KI270735v1_random//g' | sed 's/_KI270736v1_random//g' | sed 's/_KI270737v1_random//g' | sed 's/_GL000221v1_random//g' | sed 's/_GL000008v2_random//g' | sed 's/_KI270717v1_random//g' | sed 's/_KI270718v1_random//g' | sed 's/_KI270719v1_random//g' | sed 's/_KI270720v1_random//g' |sed 's/_GL383545v1_alt//g'| sed 's/_GL383546v1_alt//g'| sed 's/_KI270824v1_alt//g'| sed 's/_KI270825v1_alt//g'| sed 's/_JH159136v1_alt//g'| sed 's/_JH159137v1_alt//g'| sed 's/_KI270829v1_alt//g'| sed 's/_KI270830v1_alt//g'| sed 's/_KI270831v1_alt//g'| sed 's/_KI270832v1_alt//g'| sed 's/_KI270902v1_alt//g'| sed 's/_KI270903v1_alt//g'| sed 's/_KI270927v1_alt//g'| sed 's/_GL383550v2_alt//g'| sed 's/_GL383553v2_alt//g'| sed 's/_GL877875v1_alt//g'| sed 's/_GL877876v1_alt//g'| sed 's/_KI270833v1_alt//g'| sed 's/_KI270834v1_alt//g'| sed 's/_KI270835v1_alt//g'| sed 's/_KI270837v1_alt//g'| sed 's/_KI270904v1_alt//g'| sed 's/_KI270838v1_alt//g'| sed 's/_KI270839v1_alt//g'| sed 's/_KI270840v1_alt//g'| sed 's/_KI270842v1_alt//g'| sed 's/_KI270844v1_alt//g'| sed 's/_KI270845v1_alt//g'| sed 's/_KI270846v1_alt//g'| sed 's/_KI270847v1_alt//g'| sed 's/_GL383554v1_alt//g'| sed 's/_GL383555v2_alt//g'| sed 's/_KI270848v1_alt//g'| sed 's/_KI270849v1_alt//g'| sed 's/_KI270850v1_alt//g'| sed 's/_KI270851v1_alt//g'| sed 's/_KI270852v1_alt//g'| sed 's/_KI270905v1_alt//g'| sed 's/_KI270906v1_alt//g'| sed 's/_GL383556v1_alt//g'| sed 's/_GL383557v1_alt//g'| sed 's/_KI270853v1_alt//g'| sed 's/_KI270854v1_alt//g'| sed 's/_KI270855v1_alt//g'| sed 's/_KI270856v1_alt//g'| sed 's/_GL000258v2_alt//g'| sed 's/_GL383563v3_alt//g'| sed 's/_GL383564v2_alt//g'| sed 's/_GL383566v1_alt//g'| sed 's/_JH159146v1_alt//g'| sed 's/_JH159147v1_alt//g'| sed 's/_JH159148v1_alt//g'| sed 's/_KI270857v1_alt//g'| sed 's/_KI270858v1_alt//g'| sed 's/_KI270859v1_alt//g'| sed 's/_KI270860v1_alt//g'| sed 's/_KI270861v1_alt//g'| sed 's/_KI270862v1_alt//g'| sed 's/_KI270907v1_alt//g'| sed 's/_KI270908v1_alt//g'| sed 's/_KI270909v1_alt//g'| sed 's/_KI270910v1_alt//g'| sed 's/_GL383567v1_alt//g'| sed 's/_GL383571v1_alt//g'| sed 's/_GL383572v1_alt//g'| sed 's/_KI270863v1_alt//g'| sed 's/_KI270911v1_alt//g'| sed 's/_GL000209v2_alt//g'| sed 's/_GL383573v1_alt//g'| sed 's/_GL383574v1_alt//g'| sed 's/_GL383575v2_alt//g'| sed 's/_GL383576v1_alt//g'| sed 's/_GL949746v1_alt//g'| sed 's/_GL949747v2_alt//g'| sed 's/_GL949748v2_alt//g'| sed 's/_GL949749v2_alt//g'| sed 's/_GL949750v2_alt//g'| sed 's/_GL949751v2_alt//g'| sed 's/_GL949752v1_alt//g'| sed 's/_GL949753v2_alt//g'| sed 's/_KI270865v1_alt//g'| sed 's/_KI270866v1_alt//g'| sed 's/_KI270867v1_alt//g'| sed 's/_KI270868v1_alt//g'| sed 's/_KI270882v1_alt//g'| sed 's/_KI270883v1_alt//g'| sed 's/_KI270885v1_alt//g'| sed 's/_KI270886v1_alt//g'| sed 's/_KI270887v1_alt//g'| sed 's/_KI270890v1_alt//g'| sed 's/_KI270891v1_alt//g'| sed 's/_KI270916v1_alt//g'| sed 's/_KI270917v1_alt//g'| sed 's/_KI270918v1_alt//g'| sed 's/_KI270920v1_alt//g'| sed 's/_KI270921v1_alt//g'| sed 's/_KI270922v1_alt//g'| sed 's/_KI270923v1_alt//g'| sed 's/_KI270929v1_alt//g'| sed 's/_KI270930v1_alt//g'| sed 's/_KI270938v1_alt//g'| sed 's/_GL383518v1_alt//g'| sed 's/_GL383519v1_alt//g'| sed 's/_GL383520v2_alt//g'| sed 's/_KI270759v1_alt//g'| sed 's/_KI270760v1_alt//g'| sed 's/_KI270761v1_alt//g'| sed 's/_KI270762v1_alt//g'| sed 's/_KI270763v1_alt//g'| sed 's/_KI270766v1_alt//g'| sed 's/_KI270892v1_alt//g'| sed 's/_KI270869v1_alt//g'| sed 's/_KI270870v1_alt//g'| sed 's/_KI270871v1_alt//g'| sed 's/_GL383579v2_alt//g'| sed 's/_GL383580v2_alt//g'| sed 's/_GL383581v2_alt//g'| sed 's/_KI270872v1_alt//g'| sed 's/_KI270873v1_alt//g'| sed 's/_GL383582v2_alt//g'| sed 's/_GL383583v2_alt//g'| sed 's/_KB663609v1_alt//g'| sed 's/_KI270875v1_alt//g'| sed 's/_KI270876v1_alt//g'| sed 's/_KI270877v1_alt//g'| sed 's/_KI270878v1_alt//g'| sed 's/_KI270879v1_alt//g'| sed 's/_KI270928v1_alt//g'| sed 's/_GL383521v1_alt//g'| sed 's/_GL383522v1_alt//g'| sed 's/_GL582966v2_alt//g'| sed 's/_KI270767v1_alt//g'| sed 's/_KI270768v1_alt//g'| sed 's/_KI270769v1_alt//g'| sed 's/_KI270772v1_alt//g'| sed 's/_KI270773v1_alt//g'| sed 's/_KI270774v1_alt//g'| sed 's/_KI270776v1_alt//g'| sed 's/_KI270893v1_alt//g'| sed 's/_KI270894v1_alt//g'| sed 's/_GL383526v1_alt//g'| sed 's/_JH636055v2_alt//g'| sed 's/_KI270777v1_alt//g'| sed 's/_KI270779v1_alt//g'| sed 's/_KI270780v1_alt//g'| sed 's/_KI270781v1_alt//g'| sed 's/_KI270782v1_alt//g'| sed 's/_KI270784v1_alt//g'| sed 's/_KI270895v1_alt//g'| sed 's/_KI270924v1_alt//g'| sed 's/_KI270934v1_alt//g'| sed 's/_KI270935v1_alt//g'| sed 's/_KI270936v1_alt//g'| sed 's/_KI270937v1_alt//g'| sed 's/_GL000257v2_alt//g'| sed 's/_GL383527v1_alt//g'| sed 's/_KI270789v1_alt//g'| sed 's/_KI270790v1_alt//g'| sed 's/_KI270896v1_alt//g'| sed 's/_KI270925v1_alt//g'| sed 's/_GL339449v2_alt//g'| sed 's/_KI270791v1_alt//g'| sed 's/_KI270792v1_alt//g'| sed 's/_KI270793v1_alt//g'| sed 's/_KI270794v1_alt//g'| sed 's/_KI270795v1_alt//g'| sed 's/_KI270897v1_alt//g'| sed 's/_KI270898v1_alt//g'| sed 's/_GL000250v2_alt//g'| sed 's/_GL000251v2_alt//g'| sed 's/_GL000252v2_alt//g'| sed 's/_GL000253v2_alt//g'| sed 's/_GL000254v2_alt//g'| sed 's/_GL000255v2_alt//g'| sed 's/_GL000256v2_alt//g'| sed 's/_GL383533v1_alt//g'| sed 's/_KI270758v1_alt//g'| sed 's/_KI270797v1_alt//g'| sed 's/_KI270798v1_alt//g'| sed 's/_KI270800v1_alt//g'| sed 's/_KI270801v1_alt//g'| sed 's/_GL383534v2_alt//g'| sed 's/_KI270803v1_alt//g'| sed 's/_KI270804v1_alt//g'| sed 's/_KI270805v1_alt//g'| sed 's/_KI270806v1_alt//g'| sed 's/_KI270808v1_alt//g'| sed 's/_KI270809v1_alt//g'| sed 's/_KI270899v1_alt//g'| sed 's/_KI270811v1_alt//g'| sed 's/_KI270812v1_alt//g'| sed 's/_KI270813v1_alt//g'| sed 's/_KI270815v1_alt//g'| sed 's/_KI270816v1_alt//g'| sed 's/_KI270817v1_alt//g'| sed 's/_KI270818v1_alt//g'| sed 's/_KI270819v1_alt//g'| sed 's/_KI270820v1_alt//g'| sed 's/_KI270821v1_alt//g'| sed 's/_KI270822v1_alt//g'| sed 's/_KI270900v1_alt//g'| sed 's/_KI270926v1_alt//g'| sed 's/_GL383540v1_alt//g'| sed 's/_GL383541v1_alt//g'| sed 's/_GL383542v1_alt//g'| sed 's/_KI270823v1_alt//g'| sed 's/_KI270881v1_alt//g'| sed 's/_KI270913v1_alt//g' | gzip > PE18.txt.gz
$ zcat PE20.small.seqz.gz | sed 's/_KI270721v1_random//g'| sed 's/_GL000009v2_random//g'| sed 's/_GL000194v1_random//g'| sed 's/_GL000225v1_random//g'| sed 's/_KI270722v1_random//g'| sed 's/_KI270723v1_random//g'| sed 's/_KI270724v1_random//g'| sed 's/_KI270725v1_random//g'| sed 's/_KI270726v1_random//g'| sed 's/_KI270727v1_random//g'| sed 's/_KI270728v1_random//g'| sed 's/_GL000205v2_random//g'| sed 's/_KI270729v1_random//g'| sed 's/_KI270706v1_random//g'| sed 's/_KI270708v1_random//g'| sed 's/_KI270709v1_random//g'| sed 's/_KI270711v1_random//g'| sed 's/_KI270712v1_random//g'| sed 's/_KI270713v1_random//g'| sed 's/_KI270714v1_random//g'| sed 's/_KI270731v1_random//g'| sed 's/_KI270732v1_random//g'| sed 's/_KI270733v1_random//g'| sed 's/_KI270734v1_random//g'| sed 's/_KI270735v1_random//g'| sed 's/_KI270736v1_random//g'| sed 's/_KI270737v1_random//g'| sed 's/_GL000221v1_random//g'| sed 's/_GL000008v2_random//g'| sed 's/_KI270717v1_random//g'| sed 's/_KI270718v1_random//g'| sed 's/_KI270719v1_random//g'| sed 's/_KI270720v1_random//g'| sed 's/_GL383545v1_alt//g' | sed 's/_GL383546v1_alt//g' | sed 's/_KI270824v1_alt//g' | sed 's/_KI270825v1_alt//g' | sed 's/_JH159136v1_alt//g' | sed 's/_JH159137v1_alt//g' | sed 's/_KI270829v1_alt//g' | sed 's/_KI270830v1_alt//g' | sed 's/_KI270831v1_alt//g' | sed 's/_KI270832v1_alt//g' | sed 's/_KI270902v1_alt//g' | sed 's/_KI270903v1_alt//g' | sed 's/_KI270927v1_alt//g' | sed 's/_GL383550v2_alt//g' | sed 's/_GL383553v2_alt//g' | sed 's/_GL877875v1_alt//g' | sed 's/_GL877876v1_alt//g' | sed 's/_KI270833v1_alt//g' | sed 's/_KI270834v1_alt//g' | sed 's/_KI270835v1_alt//g' | sed 's/_KI270837v1_alt//g' | sed 's/_KI270904v1_alt//g' | sed 's/_KI270838v1_alt//g' | sed 's/_KI270839v1_alt//g' | sed 's/_KI270840v1_alt//g' | sed 's/_KI270842v1_alt//g' | sed 's/_KI270844v1_alt//g' | sed 's/_KI270845v1_alt//g' | sed 's/_KI270846v1_alt//g' | sed 's/_KI270847v1_alt//g' | sed 's/_GL383554v1_alt//g' | sed 's/_GL383555v2_alt//g' | sed 's/_KI270848v1_alt//g' | sed 's/_KI270849v1_alt//g' | sed 's/_KI270850v1_alt//g' | sed 's/_KI270851v1_alt//g' | sed 's/_KI270852v1_alt//g' | sed 's/_KI270905v1_alt//g' | sed 's/_KI270906v1_alt//g' | sed 's/_GL383556v1_alt//g' | sed 's/_GL383557v1_alt//g' | sed 's/_KI270853v1_alt//g' | sed 's/_KI270854v1_alt//g' | sed 's/_KI270855v1_alt//g' | sed 's/_KI270856v1_alt//g' | sed 's/_GL000258v2_alt//g' | sed 's/_GL383563v3_alt//g' | sed 's/_GL383564v2_alt//g' | sed 's/_GL383566v1_alt//g' | sed 's/_JH159146v1_alt//g' | sed 's/_JH159147v1_alt//g' | sed 's/_JH159148v1_alt//g' | sed 's/_KI270857v1_alt//g' | sed 's/_KI270858v1_alt//g' | sed 's/_KI270859v1_alt//g' | sed 's/_KI270860v1_alt//g' | sed 's/_KI270861v1_alt//g' | sed 's/_KI270862v1_alt//g' | sed 's/_KI270907v1_alt//g' | sed 's/_KI270908v1_alt//g' | sed 's/_KI270909v1_alt//g' | sed 's/_KI270910v1_alt//g' | sed 's/_GL383567v1_alt//g' | sed 's/_GL383571v1_alt//g' | sed 's/_GL383572v1_alt//g' | sed 's/_KI270863v1_alt//g' | sed 's/_KI270911v1_alt//g' | sed 's/_GL000209v2_alt//g' | sed 's/_GL383573v1_alt//g' | sed 's/_GL383574v1_alt//g' | sed 's/_GL383575v2_alt//g' | sed 's/_GL383576v1_alt//g' | sed 's/_GL949746v1_alt//g' | sed 's/_GL949747v2_alt//g' | sed 's/_GL949748v2_alt//g' | sed 's/_GL949749v2_alt//g' | sed 's/_GL949750v2_alt//g' | sed 's/_GL949751v2_alt//g' | sed 's/_GL949752v1_alt//g' | sed 's/_GL949753v2_alt//g' | sed 's/_KI270865v1_alt//g' | sed 's/_KI270866v1_alt//g' | sed 's/_KI270867v1_alt//g' | sed 's/_KI270868v1_alt//g' | sed 's/_KI270882v1_alt//g' | sed 's/_KI270884v1_alt//g' | sed 's/_KI270885v1_alt//g' | sed 's/_KI270886v1_alt//g' | sed 's/_KI270887v1_alt//g' | sed 's/_KI270890v1_alt//g' | sed 's/_KI270916v1_alt//g' | sed 's/_KI270917v1_alt//g' | sed 's/_KI270918v1_alt//g' | sed 's/_KI270920v1_alt//g' | sed 's/_KI270921v1_alt//g' | sed 's/_KI270922v1_alt//g' | sed 's/_KI270923v1_alt//g' | sed 's/_KI270929v1_alt//g' | sed 's/_KI270931v1_alt//g' | sed 's/_KI270938v1_alt//g' | sed 's/_GL383518v1_alt//g' | sed 's/_GL383519v1_alt//g' | sed 's/_GL383520v2_alt//g' | sed 's/_KI270759v1_alt//g' | sed 's/_KI270760v1_alt//g' | sed 's/_KI270761v1_alt//g' | sed 's/_KI270762v1_alt//g' | sed 's/_KI270763v1_alt//g' | sed 's/_KI270766v1_alt//g' | sed 's/_KI270892v1_alt//g' | sed 's/_KI270869v1_alt//g' | sed 's/_KI270870v1_alt//g' | sed 's/_KI270871v1_alt//g' | sed 's/_GL383579v2_alt//g' | sed 's/_GL383580v2_alt//g' | sed 's/_GL383581v2_alt//g' | sed 's/_KI270872v1_alt//g' | sed 's/_KI270873v1_alt//g' | sed 's/_GL383582v2_alt//g' | sed 's/_GL383583v2_alt//g' | sed 's/_KB663609v1_alt//g' | sed 's/_KI270875v1_alt//g' | sed 's/_KI270876v1_alt//g' | sed 's/_KI270877v1_alt//g' | sed 's/_KI270878v1_alt//g' | sed 's/_KI270879v1_alt//g' | sed 's/_KI270928v1_alt//g' | sed 's/_GL383521v1_alt//g' | sed 's/_GL383522v1_alt//g' | sed 's/_GL582966v2_alt//g' | sed 's/_KI270767v1_alt//g' | sed 's/_KI270768v1_alt//g' | sed 's/_KI270769v1_alt//g' | sed 's/_KI270772v1_alt//g' | sed 's/_KI270773v1_alt//g' | sed 's/_KI270774v1_alt//g' | sed 's/_KI270775v1_alt//g' | sed 's/_KI270776v1_alt//g' | sed 's/_KI270893v1_alt//g' | sed 's/_KI270894v1_alt//g' | sed 's/_GL383526v1_alt//g' | sed 's/_JH636055v2_alt//g' | sed 's/_KI270777v1_alt//g' | sed 's/_KI270779v1_alt//g' | sed 's/_KI270780v1_alt//g' | sed 's/_KI270781v1_alt//g' | sed 's/_KI270782v1_alt//g' | sed 's/_KI270784v1_alt//g' | sed 's/_KI270895v1_alt//g' | sed 's/_KI270924v1_alt//g' | sed 's/_KI270934v1_alt//g' | sed 's/_KI270935v1_alt//g' | sed 's/_KI270936v1_alt//g' | sed 's/_KI270937v1_alt//g' | sed 's/_GL000257v2_alt//g' | sed 's/_GL383527v1_alt//g' | sed 's/_KI270789v1_alt//g' | sed 's/_KI270790v1_alt//g' | sed 's/_KI270896v1_alt//g' | sed 's/_KI270925v1_alt//g' | sed 's/_GL339449v2_alt//g' | sed 's/_GL949742v1_alt//g' | sed 's/_KI270791v1_alt//g' | sed 's/_KI270792v1_alt//g' | sed 's/_KI270793v1_alt//g' | sed 's/_KI270794v1_alt//g' | sed 's/_KI270795v1_alt//g' | sed 's/_KI270897v1_alt//g' | sed 's/_KI270898v1_alt//g' | sed 's/_GL000250v2_alt//g' | sed 's/_GL000251v2_alt//g' | sed 's/_GL000252v2_alt//g' | sed 's/_GL000253v2_alt//g' | sed 's/_GL000254v2_alt//g' | sed 's/_GL000255v2_alt//g' | sed 's/_GL000256v2_alt//g' | sed 's/_GL383533v1_alt//g' | sed 's/_KI270758v1_alt//g' | sed 's/_KI270797v1_alt//g' | sed 's/_KI270798v1_alt//g' | sed 's/_KI270800v1_alt//g' | sed 's/_KI270801v1_alt//g' | sed 's/_GL383534v2_alt//g' | sed 's/_KI270803v1_alt//g' | sed 's/_KI270804v1_alt//g' | sed 's/_KI270805v1_alt//g' | sed 's/_KI270806v1_alt//g' | sed 's/_KI270808v1_alt//g' | sed 's/_KI270809v1_alt//g' | sed 's/_KI270899v1_alt//g' | sed 's/_KI270811v1_alt//g' | sed 's/_KI270812v1_alt//g' | sed 's/_KI270813v1_alt//g' | sed 's/_KI270815v1_alt//g' | sed 's/_KI270816v1_alt//g' | sed 's/_KI270817v1_alt//g' | sed 's/_KI270818v1_alt//g' | sed 's/_KI270819v1_alt//g' | sed 's/_KI270820v1_alt//g' | sed 's/_KI270821v1_alt//g' | sed 's/_KI270822v1_alt//g' | sed 's/_KI270900v1_alt//g' | sed 's/_KI270926v1_alt//g' | sed 's/_GL383540v1_alt//g' | sed 's/_GL383541v1_alt//g' | sed 's/_GL383542v1_alt//g' | sed 's/_KI270823v1_alt//g' | sed 's/_KI270880v1_alt//g' | sed 's/_KI270881v1_alt//g' | sed 's/_KI270913v1_alt//g' | gzip > PE20.txt.gz
$ zcat PE24.small.seqz.gz | sed 's/_KI270721v1_random//g' | sed 's/_GL000009v2_random//g' | sed 's/_GL000194v1_random//g' | sed 's/_GL000225v1_random//g' | sed 's/_KI270722v1_random//g' | sed 's/_KI270723v1_random//g' | sed 's/_KI270724v1_random//g' | sed 's/_KI270725v1_random//g' | sed 's/_KI270726v1_random//g' | sed 's/_KI270727v1_random//g' | sed 's/_KI270728v1_random//g' | sed 's/_GL000205v2_random//g' | sed 's/_KI270729v1_random//g' | sed 's/_KI270706v1_random//g' | sed 's/_KI270708v1_random//g' | sed 's/_KI270709v1_random//g' | sed 's/_KI270711v1_random//g' | sed 's/_KI270712v1_random//g' | sed 's/_KI270713v1_random//g' | sed 's/_KI270714v1_random//g' | sed 's/_KI270731v1_random//g' | sed 's/_KI270732v1_random//g' | sed 's/_KI270733v1_random//g' | sed 's/_KI270734v1_random//g' | sed 's/_KI270735v1_random//g' | sed 's/_KI270736v1_random//g' | sed 's/_KI270737v1_random//g' | sed 's/_GL000221v1_random//g' | sed 's/_GL000008v2_random//g' | sed 's/_KI270717v1_random//g' | sed 's/_KI270718v1_random//g' | sed 's/_KI270719v1_random//g' | sed 's/_KI270720v1_random//g' | sed 's/_GL383545v1_alt//g' | sed 's/_GL383546v1_alt//g' | sed 's/_KI270824v1_alt//g' | sed 's/_KI270825v1_alt//g' | sed 's/_JH159136v1_alt//g' | sed 's/_JH159137v1_alt//g' | sed 's/_KI270829v1_alt//g' | sed 's/_KI270830v1_alt//g' | sed 's/_KI270831v1_alt//g' | sed 's/_KI270832v1_alt//g' | sed 's/_KI270902v1_alt//g' | sed 's/_KI270903v1_alt//g' | sed 's/_KI270927v1_alt//g' | sed 's/_GL383550v2_alt//g' | sed 's/_GL383553v2_alt//g' | sed 's/_GL877875v1_alt//g' | sed 's/_GL877876v1_alt//g' | sed 's/_KI270833v1_alt//g' | sed 's/_KI270834v1_alt//g' | sed 's/_KI270835v1_alt//g' | sed 's/_KI270837v1_alt//g' | sed 's/_KI270904v1_alt//g' | sed 's/_KI270838v1_alt//g' | sed 's/_KI270839v1_alt//g' | sed 's/_KI270840v1_alt//g' | sed 's/_KI270842v1_alt//g' | sed 's/_KI270844v1_alt//g' | sed 's/_KI270845v1_alt//g' | sed 's/_KI270846v1_alt//g' | sed 's/_KI270847v1_alt//g' | sed 's/_GL383554v1_alt//g' | sed 's/_GL383555v2_alt//g' | sed 's/_KI270848v1_alt//g' | sed 's/_KI270849v1_alt//g' | sed 's/_KI270850v1_alt//g' | sed 's/_KI270851v1_alt//g' | sed 's/_KI270852v1_alt//g' | sed 's/_KI270905v1_alt//g' | sed 's/_KI270906v1_alt//g' | sed 's/_GL383556v1_alt//g' | sed 's/_GL383557v1_alt//g' | sed 's/_KI270853v1_alt//g' | sed 's/_KI270854v1_alt//g' | sed 's/_KI270855v1_alt//g' | sed 's/_KI270856v1_alt//g' | sed 's/_GL000258v2_alt//g' | sed 's/_GL383563v3_alt//g' | sed 's/_GL383564v2_alt//g' | sed 's/_GL383566v1_alt//g' | sed 's/_JH159146v1_alt//g' | sed 's/_JH159147v1_alt//g' | sed 's/_JH159148v1_alt//g' | sed 's/_KI270857v1_alt//g' | sed 's/_KI270858v1_alt//g' | sed 's/_KI270859v1_alt//g' | sed 's/_KI270860v1_alt//g' | sed 's/_KI270861v1_alt//g' | sed 's/_KI270862v1_alt//g' | sed 's/_KI270907v1_alt//g' | sed 's/_KI270908v1_alt//g' | sed 's/_KI270909v1_alt//g' | sed 's/_KI270910v1_alt//g' | sed 's/_GL383567v1_alt//g' | sed 's/_GL383571v1_alt//g' | sed 's/_GL383572v1_alt//g' | sed 's/_KI270863v1_alt//g' | sed 's/_KI270911v1_alt//g' | sed 's/_GL000209v2_alt//g' | sed 's/_GL383573v1_alt//g' | sed 's/_GL383574v1_alt//g' | sed 's/_GL383575v2_alt//g' | sed 's/_GL383576v1_alt//g' | sed 's/_GL949746v1_alt//g' | sed 's/_GL949747v2_alt//g' | sed 's/_GL949748v2_alt//g' | sed 's/_GL949749v2_alt//g' | sed 's/_GL949750v2_alt//g' | sed 's/_GL949751v2_alt//g' | sed 's/_GL949752v1_alt//g' | sed 's/_GL949753v2_alt//g' | sed 's/_KI270865v1_alt//g' | sed 's/_KI270866v1_alt//g' | sed 's/_KI270867v1_alt//g' | sed 's/_KI270868v1_alt//g' | sed 's/_KI270882v1_alt//g' | sed 's/_KI270883v1_alt//g' | sed 's/_KI270885v1_alt//g' | sed 's/_KI270886v1_alt//g' | sed 's/_KI270887v1_alt//g' | sed 's/_KI270890v1_alt//g' | sed 's/_KI270916v1_alt//g' | sed 's/_KI270917v1_alt//g' | sed 's/_KI270918v1_alt//g' | sed 's/_KI270920v1_alt//g' | sed 's/_KI270921v1_alt//g' | sed 's/_KI270922v1_alt//g' | sed 's/_KI270923v1_alt//g' | sed 's/_KI270929v1_alt//g' | sed 's/_KI270930v1_alt//g' | sed 's/_KI270938v1_alt//g' | sed 's/_GL383518v1_alt//g' | sed 's/_GL383519v1_alt//g' | sed 's/_GL383520v2_alt//g' | sed 's/_KI270759v1_alt//g' | sed 's/_KI270760v1_alt//g' | sed 's/_KI270761v1_alt//g' | sed 's/_KI270762v1_alt//g' | sed 's/_KI270763v1_alt//g' | sed 's/_KI270766v1_alt//g' | sed 's/_KI270892v1_alt//g' | sed 's/_KI270869v1_alt//g' | sed 's/_KI270870v1_alt//g' | sed 's/_KI270871v1_alt//g' | sed 's/_GL383579v2_alt//g' | sed 's/_GL383580v2_alt//g' | sed 's/_GL383581v2_alt//g' | sed 's/_KI270872v1_alt//g' | sed 's/_KI270873v1_alt//g' | sed 's/_GL383582v2_alt//g' | sed 's/_GL383583v2_alt//g' | sed 's/_KB663609v1_alt//g' | sed 's/_KI270875v1_alt//g' | sed 's/_KI270876v1_alt//g' | sed 's/_KI270877v1_alt//g' | sed 's/_KI270878v1_alt//g' | sed 's/_KI270879v1_alt//g' | sed 's/_KI270928v1_alt//g' | sed 's/_GL383521v1_alt//g' | sed 's/_GL383522v1_alt//g' | sed 's/_GL582966v2_alt//g' | sed 's/_KI270767v1_alt//g' | sed 's/_KI270768v1_alt//g' | sed 's/_KI270769v1_alt//g' | sed 's/_KI270772v1_alt//g' | sed 's/_KI270773v1_alt//g' | sed 's/_KI270774v1_alt//g' | sed 's/_KI270775v1_alt//g' | sed 's/_KI270776v1_alt//g' | sed 's/_KI270893v1_alt//g' | sed 's/_KI270894v1_alt//g' | sed 's/_GL383526v1_alt//g' | sed 's/_JH636055v2_alt//g' | sed 's/_KI270777v1_alt//g' | sed 's/_KI270779v1_alt//g' | sed 's/_KI270780v1_alt//g' | sed 's/_KI270781v1_alt//g' | sed 's/_KI270782v1_alt//g' | sed 's/_KI270784v1_alt//g' | sed 's/_KI270895v1_alt//g' | sed 's/_KI270924v1_alt//g' | sed 's/_KI270934v1_alt//g' | sed 's/_KI270935v1_alt//g' | sed 's/_KI270936v1_alt//g' | sed 's/_KI270937v1_alt//g' | sed 's/_GL000257v2_alt//g' | sed 's/_GL383527v1_alt//g' | sed 's/_KI270789v1_alt//g' | sed 's/_KI270790v1_alt//g' | sed 's/_KI270896v1_alt//g' | sed 's/_KI270925v1_alt//g' | sed 's/_GL339449v2_alt//g' | sed 's/_KI270791v1_alt//g' | sed 's/_KI270792v1_alt//g' | sed 's/_KI270793v1_alt//g' | sed 's/_KI270794v1_alt//g' | sed 's/_KI270795v1_alt//g' | sed 's/_KI270897v1_alt//g' | sed 's/_KI270898v1_alt//g' | sed 's/_GL000250v2_alt//g' | sed 's/_GL000251v2_alt//g' | sed 's/_GL000252v2_alt//g' | sed 's/_GL000253v2_alt//g' | sed 's/_GL000254v2_alt//g' | sed 's/_GL000255v2_alt//g' | sed 's/_GL000256v2_alt//g' | sed 's/_GL383533v1_alt//g' | sed 's/_KI270758v1_alt//g' | sed 's/_KI270797v1_alt//g' | sed 's/_KI270798v1_alt//g' | sed 's/_KI270800v1_alt//g' | sed 's/_KI270801v1_alt//g' | sed 's/_GL383534v2_alt//g' | sed 's/_KI270803v1_alt//g' | sed 's/_KI270804v1_alt//g' | sed 's/_KI270805v1_alt//g' | sed 's/_KI270806v1_alt//g' | sed 's/_KI270808v1_alt//g' | sed 's/_KI270809v1_alt//g' | sed 's/_KI270899v1_alt//g' | sed 's/_KI270811v1_alt//g' | sed 's/_KI270812v1_alt//g' | sed 's/_KI270813v1_alt//g' | sed 's/_KI270815v1_alt//g' | sed 's/_KI270816v1_alt//g' | sed 's/_KI270817v1_alt//g' | sed 's/_KI270818v1_alt//g' | sed 's/_KI270819v1_alt//g' | sed 's/_KI270820v1_alt//g' | sed 's/_KI270821v1_alt//g' | sed 's/_KI270822v1_alt//g' | sed 's/_KI270900v1_alt//g' | sed 's/_KI270926v1_alt//g' | sed 's/_GL383540v1_alt//g' | sed 's/_GL383541v1_alt//g' | sed 's/_GL383542v1_alt//g' | sed 's/_KI270823v1_alt//g' | sed 's/_KI270881v1_alt//g' | gzip > PE24.txt.gz

test <- sequenza.extract_hg38("/storage2/Project/CSC/WES/04_CNV/Sequenza/w_50/PE17.txt.gz")

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
PyClone build_mutations_file TSV_FILE