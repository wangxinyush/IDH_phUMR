#annoation of partial hyper, full hyper, flanking umr

blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH GBM
red = rgb(255, 45, 41, max=255) #WT GBM
purple = rgb(153, 0, 153, max=255) #
green = rgb(12, 99, 52, max=255) #CGI

cgi_blue = rgb(57, 127, 185, max=255) #CGI
shore_purple = rgb(152, 79, 159, max=255) #shore
shelf_green = rgb(80, 175, 73, max=255) #shelf
sea_gray = rgb(137, 137, 137, max=255) #sea

#Data
partial_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\sort_partial_hyper.bed",
                            header=F, sep="\t", fill=T, check.names=F, quote="")
full_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\sort_full_hyper.bed",
                         header=F, sep="\t", fill=T, check.names=F, quote="")
flanking_umr <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\sort_partial_adj_umr.bed",
                           header=F, sep="\t", fill=T, check.names=F, quote="")

#1) basic stat - length
#par(mar=c(5.1,4.1,4.1,2.1))
partial_hyper_length <- partial_hyper$V3 - partial_hyper$V2 + 1
full_hyper_length <- full_hyper$V3 - full_hyper$V2 + 1
flanking_umr_length <- flanking_umr$V3 - flanking_umr$V2 + 1

boxplot(full_hyper_length, partial_hyper_length, flanking_umr_length,
        col=c(orange, red, cgi_blue), las=1, bty="l", outline = F,
        cex.axis = 1.5, cex.lab=1.5, cex.main = 2, 
        names=c("Fully", "Partially", "UMR"))

plot(density( partial_hyper_length ) ,
     xlab = "Length", ylab = "Density",
     col=red, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
lines(density(full_hyper_length), ylab="", xlab="", col=orange, lwd=2, cex.lab=1.3, las=1)
lines(density(flanking_umr_length), ylab="", xlab="", col=cgi_blue, lwd=2, cex.lab=1.3, las=1)

dev.print(pdf, file="Rplot/parital_full_hyper_umr_length.pdf", width=4, height=6)
dev.print(pdf, file="Rplot/parital_full_hyper_umr_length_density.pdf", width=4, height=4)
#legend(2000,0.001,text.font=2, c("partial", "full"), lwd=c(2,2), lty=c(1,1), col=c(red,"gray")) 

#2) genomic distribution
library("ChIPseeker")
library("GenomicRanges")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

partial_hyper_bed <- GRanges(seqnames=partial_hyper[,1], ranges=IRanges(start=partial_hyper[,2], end=partial_hyper[,3]))
full_hyper_bed <- GRanges(seqnames=full_hyper[,1], ranges=IRanges(start=full_hyper[,2], end=full_hyper[,3]))
flanking_umr_bed <- GRanges(seqnames=flanking_umr[,1], ranges=IRanges(start=flanking_umr[,2], end=flanking_umr[,3]))

full_Hyper_anno <- annotatePeak(full_hyper_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
full_Hyper_anno_table <- as.data.frame(full_Hyper_anno)
partial_hyper_anno <- annotatePeak(partial_hyper_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
partial_hyper_anno_table <- as.data.frame(partial_hyper_anno)
flanking_umr_anno <- annotatePeak(flanking_umr_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
flanking_umr_anno_table <- as.data.frame(flanking_umr_anno)

anno_list <- list(full_Hyper_anno_table$annotation, partial_hyper_anno_table$annotation, flanking_umr_anno_table$annotation)
anno_names <- c("full_hyper", "partial_hyper", "flanking_umr")
anno_color <- c(orange, red, cgi_blue)

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
        ylab="Percent", ylim=c(0, 1))
dev.print(pdf, file="Rplot/partial_full_hyper_umr_annotation.pdf", width=4, height=5)
par(mar=c(5.1,4.1,4.1,2.1))

#3) CGI overlap
partial_CGI_anno <- read.table("CGI_anno/partial_hyper_CGI_anno.txt", header=F, sep="\t", fill=T, check.names=F, quote="")
full_CGI_anno <- read.table("CGI_anno/full_hyper_CGI_anno.txt", header=F, sep="\t", fill=T, check.names=F, quote="")
partial_umr_CGI_anno <- read.table("CGI_anno/partial_adj_umr_CGI_anno.txt", header=F, sep="\t", fill=T, check.names=F, quote="")

cgi_blue = rgb(57, 127, 185, max=255) #CGI
shore_purple = rgb(152, 79, 159, max=255) #shore
shelf_green = rgb(80, 175, 73, max=255) #shelf
sea_gray = rgb(137, 137, 137, max=255) #sea

CGI_percent <- table(partial_CGI_anno$V4) / nrow(partial_CGI_anno)
CGI_percent <- table(full_CGI_anno$V4) / nrow(full_CGI_anno)
CGI_percent <- table(partial_umr_CGI_anno$V4) / nrow(partial_umr_CGI_anno)

names(CGI_percent) <- c(paste("CGI shelf",paste(round( CGI_percent[1]*100, 2), "%", sep=""), sep=" "),
                        paste("CGI shore",paste(round( CGI_percent[2]*100, 2), "%", sep=""), sep=" "),
                        paste("CpG islands",paste(round( CGI_percent[3]*100, 2), "%", sep=""), sep=" "),
                        paste("Open sea",paste(round( CGI_percent[4]*100, 2), "%", sep=""), sep=" ")
)

pie(CGI_percent, col = c(shelf_green, shore_purple, cgi_blue, sea_gray), border="white")

dev.print(pdf, file="Rplot/partial_hyper_CGI.pdf", width=8, height=8)
dev.print(pdf, file="Rplot/full_hyper_CGI.pdf", width=8, height=8)
dev.print(pdf, file="Rplot/parital_umr_CGI.pdf", width=8, height=8)

#4) Methylation changes
partial_methy_matrix <- read.table("methy_matrix\\sort_partial_hyper_methy_matrix.txt", 
                                   header=T, sep="\t", fill=T, check.names=F, quote="")
full_methy_matrix <- read.table("methy_matrix\\sort_full_hyper_methy_matrix.txt", 
                                header=T, sep="\t", fill=T, check.names=F, quote="")
umr_methy_matrix <- read.table("methy_matrix\\sort_partial_umr_methy_matrix.txt", 
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

umr_methy_changes <- c()
for(i in 1:nrow(umr_methy_matrix)){
  normal_methy <- na.omit(as.numeric(umr_methy_matrix[i, normal_index]))
  IDH_methy <- na.omit(as.numeric(umr_methy_matrix[i, IDH_index]))
  WT_methy <- na.omit(as.numeric(umr_methy_matrix[i, WT_index]))
  umr_methy_changes <- c(umr_methy_changes, mean(IDH_methy) - mean(normal_methy))
}

boxplot(full_methy_changes, partial_methy_changes, umr_methy_changes,
        col=c(orange, red, cgi_blue ), las=1, bty="l", outline = F,
        cex.axis = 1.5, cex.lab=1.5, cex.main = 2, ylim=c(0,1),
        names=c("Full_Hyper", "Partially_Hyper", "Partially_adj_UMR"),
        main=t.test(partial_methy_changes, full_methy_changes)$p.value)
dev.print(pdf, file="Rplot/partial_full_hyper_umr_methy_changes.pdf", width=6, height=8)


#5) Conservation score
partial_hyper_cons <- read.table("cons\\partial_hyper_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
partial_adj_umr_cons <- read.table("cons\\partial_adj_umr_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")
full_hyper_cons <- read.table("cons\\full_hyper_cons\\region_cons.bed", header=T, sep="\t", fill=T, check.names=F, quote="")

boxplot(full_hyper_cons$mean_cons, partial_hyper_cons$mean_cons, partial_adj_umr_cons$mean_cons,
        col=c(orange, red, cgi_blue ), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 1, outline = F, ylim = c(0, 0.7),
        names=c("Full_Hyper", "Partially_Hyper", "Partially_adj_UMR"),
        main=paste(round(t.test(full_hyper_cons$mean_cons, partial_hyper_cons$mean_cons)$p.value,5),
                   t.test(partial_hyper_cons$mean_cons, partial_adj_umr_cons$mean_cons)$p.value, sep = ";") )
dev.print(pdf, file="Rplot/partial_full_hyper_umr_cons.pdf", width=6, height=8)

