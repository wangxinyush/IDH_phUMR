#features of phUMRs and fhUMRs

blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH GBM
red = rgb(255, 45, 41, max=255) #WT GBM
purple = rgb(153, 0, 153, max=255) #
green = rgb(12, 99, 52, max=255) #CGI

#Data
partial_hyper <- read.table("E:\\IDH_Hyper\\Results\\1_Identificaiton_IDH_Hyper\\202203_compared_with_CGI\\Data\\sort_partial_refumrs_hyper.bed",
                            header=F, sep="\t", fill=T, check.names=F, quote="")
full_hyper <- read.table("E:\\IDH_Hyper\\Results\\1_Identificaiton_IDH_Hyper\\202203_compared_with_CGI\\Data\\sort_full_refumrs_hyper.bed",
                            header=F, sep="\t", fill=T, check.names=F, quote="")
#1) basic stat - length
#par(mar=c(5.1,4.1,4.1,2.1))
plot(density( full_hyper$V3 - full_hyper$V2 + 1) ,
     xlab = "refUMR length", ylab = "Density",
     col="gray", lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
lines(density( partial_hyper$V3 - partial_hyper$V2 + 1), ylab="", xlab="", col=red, lwd=2, cex.lab=1.3, las=1)
dev.print(pdf, file="Rplot/parital_full_refumr_length.pdf", width=4, height=4)
#legend(2000,0.001,text.font=2, c("partial", "full"), lwd=c(2,2), lty=c(1,1), col=c(red,"gray")) 

#2) genomic distribution
library("ChIPseeker")
library("GenomicRanges")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

partial_refumr_bed <- GRanges(seqnames=partial_hyper[,1], ranges=IRanges(start=partial_hyper[,2], end=partial_hyper[,3]))
full_refumr_bed <- GRanges(seqnames=full_hyper[,1], ranges=IRanges(start=full_hyper[,2], end=full_hyper[,3]))

partial_IDH_Hyper_anno <- annotatePeak(partial_refumr_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
partial_IDH_Hyper_anno_table <- as.data.frame(partial_IDH_Hyper_anno)
full_IDH_Hyper_anno <- annotatePeak(full_refumr_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
full_IDH_Hyper_anno_table <- as.data.frame(full_IDH_Hyper_anno)

anno_list <- list(partial_IDH_Hyper_anno_table$annotation, full_IDH_Hyper_anno_table$annotation)
anno_names <- c("partial_refumr", "full_refumr")
anno_color <- c(red, "gray")

anno_per_matrix <- c()
for(i in 1:length(anno_list)){
  anno <- anno_list[[i]]
  count <- c( 
    length(grep("^Promoter", anno)),
    length(grep("^Exon", anno)) + length(grep("^5' UTR", anno)) + length(grep("^3' UTR", anno)), #Exon
    length(grep("^Intron", anno)),
    length(grep("^Distal Intergenic", anno)) + length(grep("^Downstream", anno))
  ) 
  per <- count / length(anno)
  anno_per_matrix <- rbind(anno_per_matrix, per)
}
rownames(anno_per_matrix) <- anno_names
colnames(anno_per_matrix) <- c("Promoter", "Exon", "Intron", "Intergenic")

par(mar=c(6.1,4.1,2.1,2.1))
barplot(anno_per_matrix, beside=T, col=anno_color, border=T, 
        las=2, cex.names=1.3, cex.lab=1.3, cex.axis=1.3, 
        ylab="Percent", ylim=c(0, 0.8))
dev.print(pdf, file="Rplot/partial_full_refumr_annotation.pdf", width=4, height=5)
par(mar=c(5.1,4.1,4.1,2.1))

#ggplot(data=anno_per_matrix)

#3) TSS distance
tss_dis_list <- list(partial_IDH_Hyper_anno_table$distanceToTSS, full_IDH_Hyper_anno_table$distanceToTSS)
tss_dis_names <- c("partial_refumr", "full_refumr")
tss_dis_color <- c(red, "gray")

tss_dis_per_matrix <- c()
for(i in 1:length(tss_dis_list)){
  tss_dis <- tss_dis_list[[i]]
  count <- c(
    length(which(tss_dis < -5000)),
    length(which(tss_dis >= -5000 & tss_dis < -2000)),
    length(which(tss_dis >= -2000 & tss_dis < -1000)),
    length(which(tss_dis >= -1000 & tss_dis < 0)),
    length(which(tss_dis == 0)),
    length(which(tss_dis > 0 & tss_dis <= 1000)),
    length(which(tss_dis > 1000 & tss_dis <= 2000)),
    length(which(tss_dis > 2000 & tss_dis <= 5000)),
    length(which(tss_dis > 5000))
  ) 
  per <- count / length(tss_dis)
  tss_dis_per_matrix <- rbind(tss_dis_per_matrix, per)
}
rownames(tss_dis_per_matrix) <- anno_names
colnames(tss_dis_per_matrix) <- c("<-5", "-5--2", "-2--1", "-1-0","0","0-1","1-2","2-5",">5")
barplot(tss_dis_per_matrix, beside=T, col=anno_color, border=T, las=1, cex.lab=1.3, ylab="Percent", ylim=c(0, 0.8))

#4) CGI overlap
partial_CGI_anno <- read.table("partial_refumrs_CGI_anno.txt", header=F, sep="\t", fill=T, check.names=F, quote="")
full_CGI_anno <- read.table("full_refumrs_CGI_anno.txt", header=F, sep="\t", fill=T, check.names=F, quote="")

CGI_percent <- table(partial_CGI_anno$V4) / nrow(partial_CGI_anno)
CGI_percent <- table(full_CGI_anno$V4) / nrow(full_CGI_anno)
#names(CGI_percent) : c("CGI_shelf", "CGI_shore", "CpG_islands", "Open_sea")
names(CGI_percent) <- c(paste("CGI shelf",paste(round( CGI_percent[1]*100, 2), "%", sep=""), sep=" "),
                     paste("CGI shore",paste(round( CGI_percent[2]*100, 2), "%", sep=""), sep=" "),
                     paste("CpG islands",paste(round( CGI_percent[3]*100, 2), "%", sep=""), sep=" "),
                     paste("Open sea",paste(round( CGI_percent[4]*100, 2), "%", sep=""), sep=" ")
                   )

cgi_blue = rgb(57, 127, 185, max=255) #CGI
shore_purple = rgb(152, 79, 159, max=255) #shore
shelf_green = rgb(80, 175, 73, max=255) #shelf
sea_gray = rgb(137, 137, 137, max=255) #sea

pie(CGI_percent, col = c(shelf_green, shore_purple, cgi_blue, sea_gray), border="white")

dev.print(pdf, file="Rplot/partial_refumr_CGI.pdf", width=8, height=8)
dev.print(pdf, file="Rplot/full_refumr_CGI.pdf", width=8, height=8)
#par(mfrow=c(1,1))

#5) methy changes
partial_methy_matrix <- read.table("methy_matrix\\sort_partial_hyper_methy_matrix.txt", 
                                   header=T, sep="\t", fill=T, check.names=F, quote="")
full_methy_matrix <- read.table("methy_matrix\\sort_full_hyper_methy_matrix.txt", 
                                header=T, sep="\t", fill=T, check.names=F, quote="")

sample_name <- c("chr", "start", "end", rep("Normal", 75), rep("WT", 66))
sample_name[c(83, 86, 92, 93, 97, 100, 107, 108, 111, 129, 131, 137, 139, 140, 141)] <- "IDH"
#check - right
normal_index <- c(4:78)
IDH_index <- which(sample_name == "IDH")
WT_index <- which(sample_name == "WT")

partial_methy_changes <- c()
for(i in 1:nrow(partial_methy_matrix)){
  normal_methy <- na.omit(as.numeric(partial_methy_matrix[i,normal_index]))
  IDH_methy <- na.omit(as.numeric(partial_methy_matrix[i,IDH_index]))
  WT_methy <- na.omit(as.numeric(partial_methy_matrix[i,WT_index]))
  partial_methy_changes <- c(partial_methy_changes, mean(IDH_methy) - mean(normal_methy))
}

full_methy_changes <- c()
for(i in 1:nrow(full_methy_matrix)){
  normal_methy <- na.omit(as.numeric(full_methy_matrix[i,normal_index]))
  IDH_methy <- na.omit(as.numeric(full_methy_matrix[i,IDH_index]))
  WT_methy <- na.omit(as.numeric(full_methy_matrix[i,WT_index]))
  full_methy_changes <- c(full_methy_changes, mean(IDH_methy) - mean(normal_methy))
}


boxplot(partial_methy_changes, full_methy_changes,
        col=c(red, "gray"), las=1, bty="l", outline = F,
        cex.axis = 1.5, cex.lab=1.5, cex.main = 2, ylim=c(0,1),
        names=c("Partially", "Fully"),
        main=t.test(partial_methy_changes, full_methy_changes)$p.value)
dev.print(pdf, file="Rplot/partial_full_methy_changes.pdf", width=6, height=8)

#6) conservation score
partial_refumr_cons <- read.table("cons\\partial_refumr_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
full_refumr_cons <- read.table("cons\\full_refumr_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
partial_hyper_cons <- read.table("cons\\partial_hyper_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
partial_adj_umr_cons <- read.table("cons\\partial_adj_umr_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
full_hyper_cons <- read.table("cons\\full_hyper_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")

boxplot(partial_refumr_cons$mean_cons, full_refumr_cons$mean_cons, 
        partial_hyper_cons$mean_cons, partial_adj_umr_cons$mean_cons,
        full_hyper_cons$mean_cons,
        col=c(red, "gray", shore_purple, cgi_blue, orange), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 1, outline = F, ylim = c(0, 0.7),
        names=c("Partially", "Fully", "Partially_Hyper", "Partially_adj_UMR", "Full_Hyper"),
        main=paste(round(t.test(partial_refumr_cons$mean_cons, full_refumr_cons$mean_cons)$p.value,5),
                   t.test(partial_hyper_cons$mean_cons, partial_adj_umr_cons$mean_cons)$p.value, sep = ";") )
dev.print(pdf, file="Rplot/partial_full_cons.pdf", width=6, height=8)
