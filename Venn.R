setwd("/storage2/Project/CSC/WES/03_SNV/Strelka")

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
        MuTect = c('19:8895716', '17:30563283', '17:41167816', '19:38517363', '19:55482423', '5:141151754', '16:31130933', '6_GL000255v2_alt:2611899', '7:78019354', '19:41090053'),
        MuTect2 = c('17:30563283', '2:162277587', '2:151275833', '19:37885289', '16:1979313', '5:141123569', '3:45901156', '4:4226683', '14:92083244', '2:162277594'),
        Strelka = c('20:290762', '6:129315856', '6:70859570', '19:37293857', '12:54183827', '16:46387260', 'Un_KI270750v1:92030', '6:145671385', '11:48136933', '15:22372693', 'Un_KI270467v1:2472', '20:1670956', '7:92241194', '12:70532041', 'Un_KI270442v1:223677', '10:27051454', '21:5230670', '18:11864648', '19:28791558', '13:63832928', '13:63832927', '13:63832925', '12:123321810', '22:11628394', '20:63701872', '13:21176399', '12:30983237', '15:20669797', '19:3979419', '2:27484252', '6:43784758', '14_GL000225v1_random:137242', '20:3800220', '17:42709902', 'Un_KI270744v1:107471', '14:101660331', 'Un_GL000224v1:3497', '19:28791561', '19:28791562', '9_KI270719v1_random:164235', '6:39814503', '2:167258810', '4:151722160', '4:49120192', '4:49120190', '12:132241485', '2:240871491', '22:23295013', '17:30563283', 'Un_GL000224v1:48282', '20:3800218', '4:70764506', '19:35742299', '8:64616285', '8:37965770', '4:49120189', '1:156959040', '8:68080512', '11:35301427', '19:23465926', '17_GL000205v2_random:57007', '11:99819697', '7:130310875', '10:103175705', '5:49600516', '5:13786396', '21:32334219', '2:68157216', '22:43245419', '12:119156536', '16:46391933', '22:25231937', 'Un_KI270744v1:110735', '2:209993352', '13:49691102', '11:64713239', '3:126863871', '1:154429468', 'X:101883478', '6:54380409', '7:78019354', '14:91340066', '3:16604591', '1:6234421', '11:7039639', '19:38517363', '16:34582168', '22:20405420', '22:25760998', '20:37944117', '20:30815805', '5:128536495', '5:49600535', '1:228458687', '22:35263300', '1:143232082', '10:17129922', '17:6999314', '15:22241576', '10:108299326', '1:229689691', '2:241504008', '1:143255670', '12:101369546', '15:22372723', '19:55482423', '9:21209571', '11:130911872', '9:21209574', '3:15645231', '1:29317989', '4:185408335', '5:123636251', 'X:53580794', '17:21850559', '20:30813929', 'Un_KI270467v1:2415', '14:102922666', '15:101764762', '20:31074395', '9_KI270719v1_random:163036', 'Un_KI270442v1:391855', '12:21527311', 'X:35626252', '7:152152745', '14_GL000225v1_random:58990', '19:8604984', '11:133354270', 'Un_KI270442v1:380449', '21:10414455', 'X:49174398', '1:187715296', '15:89633090', '4:49096243', '5:172670304', '19:40881764', '7:34939286', '8:144053168', '6:149541097', '8:104833514', '14_GL000225v1_random:66267', 'Un_KI270744v1:107465', '21:10416135')
        ),
    imagetype = 'png',
    filename = "PE17_M1M2S.png",
    main = "PE17 SNV call",
    main.pos = c(0.5,1),
    main.fontface = "bold",
    main.fontfamily = "serif",
    main.cex = 2,
    col = "black",
    lty = "dotted",
    lwd = 1,
    fill = c("cornflowerblue", "red",'green'),
    alpha = 0.50,
    label.col = c("white", "black", "white", "black" ,"black", "black", "white"),
    cex = 2,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkred","darkgreen"),
    cat.cex = 1.5,
    cat.fontfamily = "serif"
)

########################
library(eulerr)
library(magrittr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

setwd("/storage2/Project/CSC/WES/03_SNV/Strelka")

MuTect = c('3:143973559', '17:50201473', '19:55482423', '22:31874321', '17:50186488', '7:100953702', '17:50185889', '16:70130430', '17:50187101', '19:38517363', '10:51698416', '4:150849456', '17:50190895', '17:50186494', '10:100993347', '5:33452400', '17:50201466', '7:78019354', '17:50185904', '12:6064266', '19:41090053', '19:35861972', '4:70644248')
MuTect2 = c('13:107866037', '17:50197046', '10:119676631', '6:121447453', '21:33539220', '17:50190859', '17:50195600', '1:150151009', '4:119026703', 'X:129511850', '11:126345417', '17:50197028', '4:150849456', '17:50190895', '19:40214007', '19:15170800', '19:57859906', '10:100993347', '1:159952379', '17:50191398', '3:48433068', '1:13371759', '12:6064266', '19:57859891', '19:35861972', '17:50195603')
Strelka = c('8:143468408', )

VennDiag <- euler(c(
    "MuTect" = length(MuTect),
    "MuTect2" = length(MuTect2), 
    "Strelka" = length(Strelka), 
    "MuTect&MuTect2" = intersect(MuTect,MuTect2) %>% length, 
    "MuTect2&Strelka" = intersect(MuTect2,Strelka) %>% length, 
    "MuTect&Strelka" = intersect(MuTect,Strelka) %>% length, 
    "MuTect&MuTect2&Strelka" = MuTect %>% intersect(MuTect2) %>% intersect(Strelka) %>% length))

png("PE32_M1M2S_venn.png", width=300, height=300 )
plot1 <- plot(VennDiag, quantities = TRUE , labels = c('MuTect','MuTect2','Strelka'),fontsize = 100)
grid.arrange(grobs = list(plot1), top=textGrob("PE32", gp=gpar(fontsize=25)))
dev.off()
