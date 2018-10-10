setwd("/storage2/Project/CSC/WES/03_SNV/VarScan2")

source("https://bioconductor.org/biocLite.R")
biocLite("VennDiagram")
library(VennDiagram)
최대 5개까지 그릴 수 있음

PE17 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE17.txt", header =T, stringsAsFactors = F)
PE18 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE18.txt", header =T, stringsAsFactors = F)
PE20 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE20.txt", header =T, stringsAsFactors = F)
PE24 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE24.txt", header =T, stringsAsFactors = F)
PE32 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE32.txt", header =T, stringsAsFactors = F)
PE36 <- read.delim("/storage2/Project/CSC/WES/03_SNV/VarScan2/PE36.txt", header =T, stringsAsFactors = F)


venn.diagram(
    x = list(
        PE17 = PE17[,1],
        PE18 = PE18[,1],
        PE20 = PE20[,1],
        PE24 = PE24[,1],
        PE32 = PE32[,1],
        PE36 = PE36[,1]
        ),
    filename = "Varscan2Somatic_Patient01.png",
    col = "black",
    lty = "dotted",
    lwd = 4,
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1", "cornflowerblue", "green"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4","cornflowerblue", "green"),
    cat.cex = 2.5,
    cat.fontfamily = "serif"
)



venn.diagram(
    x = list(
        MuTect = c('5:141152011','5:141214597','7:78019354','8:30860054','13:19426473','13:19426474','13:19426481','13:44949744','15:82434534','15:82434538','16:70130430','19:36017444','19:38517363','19:41090053','19:55482423'),
        MuTect2 = c('7:100955573','11:95025680','11:95025682','12:117226688','13:19426481','13:44949744','14:68882543','19:8898976','19:36017444')
        ),
    filename = "PE24.pdf",
    main = "PE24 SNV call",
    main.pos = c(0.5,1),
    main.fontface = "bold",
    main.fontfamily = "serif",
    main.cex = 2.7,
    col = "black",
    lty = "dotted",
    lwd = 1,
    fill = c("cornflowerblue", "red"),
    alpha = 0.50,
    label.col = c("white", "black", "white"),
    cex = 2.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.pos = c(0, 0),
    cat.col = c("darkblue", "darkred"),
    cat.cex = 2,
    cat.fontfamily = "serif"
)

