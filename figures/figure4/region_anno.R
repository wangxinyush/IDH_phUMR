
setwd("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202208_oncogenes_Track\\pairs_hist_diff\\get_umr_hyper")

##1) 区域注释
# brain_glioma <- read.table("umr_200_mean0.1\\brain149_AK076_umr_hyper.bed",
#                              header=T, sep="\t", fill=T, check.names=F, quote="")
brain_glioma <- read.table("umr_200_mean0.1\\brain150_AK213_umr_hyper.bed",
                           header=T, sep="\t", fill=T, check.names=F, quote="")

length(which(brain_glioma$hyper_type == "Partially"))

partial_hyper <- brain_glioma[which(brain_glioma$hyper_type == "Partially"),]
full_hyper <- brain_glioma[which(brain_glioma$hyper_type == "Fully"),]

library("ChIPseeker")
library("GenomicRanges")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

partial_hyper_bed <- GRanges(seqnames=partial_hyper[,1], ranges=IRanges(start=partial_hyper[,5], end=partial_hyper[,6]))
full_hyper_bed <- GRanges(seqnames=full_hyper[,1], ranges=IRanges(start=full_hyper[,2], end=full_hyper[,3]))

partial_hyper_anno <- annotatePeak(partial_hyper_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
partial_hyper_anno_table <- as.data.frame(partial_hyper_anno)

full_hyper_anno <- annotatePeak(full_hyper_bed, tssRegion = c(-2000, 2000), TxDb = txdb, annoDb="org.Hs.eg.db", level="gene")
full_hyper_anno_table <- as.data.frame(full_hyper_anno)

# write.table( partial_hyper_anno_table[which(partial_hyper_anno_table$annotation != "Distal Intergenic"), ], 
#              "region_anno\\brain149_AK076_partial_hyper_gene_anno.txt", sep="\t", row.names=T, col.names=T, quote=F )
# write.table( full_hyper_anno_table[which(full_hyper_anno_table$annotation != "Distal Intergenic"), ], 
#              "region_anno\\brain149_AK076_full_hyper_gene_anno.txt", sep="\t", row.names=T, col.names=T, quote=F )

write.table( partial_hyper_anno_table[which(partial_hyper_anno_table$annotation != "Distal Intergenic"), ], 
             "region_anno\\brain150_AK213_partial_hyper_gene_anno.txt", sep="\t", row.names=T, col.names=T, quote=F )
write.table( full_hyper_anno_table[which(full_hyper_anno_table$annotation != "Distal Intergenic"), ], 
             "region_anno\\brain150_AK213_full_hyper_gene_anno.txt", sep="\t", row.names=T, col.names=T, quote=F )

##2) 基因特异性分析
#有多少是只有partially hyper, 有多少是既有partially hyper又有fully hyper的?
pair1_partial_anno <- read.table("region_anno\\brain149_AK076_partial_hyper_gene_anno.txt",
                           header=T, sep="\t", fill=T, check.names=F, quote="")
pair1_full_anno <- read.table("region_anno\\brain149_AK076_full_hyper_gene_anno.txt",
                           header=T, sep="\t", fill=T, check.names=F, quote="")

pair2_partial_anno <- read.table("region_anno\\brain150_AK213_partial_hyper_gene_anno.txt",
                                 header=T, sep="\t", fill=T, check.names=F, quote="")
pair2_full_anno <- read.table("region_anno\\brain150_AK213_full_hyper_gene_anno.txt",
                              header=T, sep="\t", fill=T, check.names=F, quote="")

#brain149_AK076, 509个重叠
length( unique(pair1_partial_anno$SYMBOL) ) #2476
length( unique(pair1_full_anno$SYMBOL) ) #3201
length( intersect( unique(pair1_partial_anno$SYMBOL),  unique(pair1_full_anno$SYMBOL) ) ) #509
#brain150_AK213, 504个重叠
length( unique(pair2_partial_anno$SYMBOL) ) #2634
length( unique(pair2_full_anno$SYMBOL) ) #3023
length( intersect( unique(pair2_partial_anno$SYMBOL),  unique(pair2_full_anno$SYMBOL) ) ) #504

#Promoter
#partial_promoter_index - partial - promoter 在pair1_partial_anno里是Promoter
#partial_body_index - partial - body, 在pair1_partial_anno里不是Promoter
#full_promoter_index - full - promoter, 在pair1_full_anno里是promoter
#full_body_index - full -body , pair1_full_anno
partial_promoter_index <- grep("Promoter", pair1_partial_anno$annotation)
partial_body_index <- grep("Promoter", pair1_partial_anno$annotation, invert = T)
full_promoter_index <- grep("Promoter", pair1_full_anno$annotation)
full_body_index <- grep("Promoter", pair1_full_anno$annotation, invert = T)

#promoter/body-partial, 没有full hyper
#promoter/body-full, 没有partial hyper
pair1_gene_list <- list(
  setdiff(unique(pair1_partial_anno$SYMBOL[partial_promoter_index]), unique(pair1_full_anno$SYMBOL) ), #promoter partial
  setdiff(unique(pair1_partial_anno$SYMBOL[partial_body_index]), unique(pair1_full_anno$SYMBOL) ), #body partial
  setdiff(unique(pair1_full_anno$SYMBOL[full_promoter_index]), unique(pair1_partial_anno$SYMBOL) ), #promoter full
  setdiff(unique(pair1_full_anno$SYMBOL[full_body_index]), unique(pair1_partial_anno$SYMBOL) ) #body full
)
names(pair1_gene_list) <- c("partial_promoter", "partial_body", "full_promoter", "full_body")

partial_promoter_index <- grep("Promoter", pair2_partial_anno$annotation)
partial_body_index <- grep("Promoter", pair2_partial_anno$annotation, invert = T)
full_promoter_index <- grep("Promoter", pair2_full_anno$annotation)
full_body_index <- grep("Promoter", pair2_full_anno$annotation, invert = T)
pair2_gene_list <- list(
  setdiff(unique(pair2_partial_anno$SYMBOL[partial_promoter_index]), unique(pair2_full_anno$SYMBOL) ), #promoter partial
  setdiff(unique(pair2_partial_anno$SYMBOL[partial_body_index]), unique(pair2_full_anno$SYMBOL) ), #body partial
  setdiff(unique(pair2_full_anno$SYMBOL[full_promoter_index]), unique(pair2_partial_anno$SYMBOL) ), #promoter full
  setdiff(unique(pair2_full_anno$SYMBOL[full_body_index]), unique(pair2_partial_anno$SYMBOL) ) #body full
)
names(pair2_gene_list) <- c("partial_promoter", "partial_body", "full_promoter", "full_body")



write_promoter_bed <- function(symbols, outfile, upstream, downstream, type){ #for bigWigAverageOverBed
  #symbols <- head(Potential_DEG_8027$symbol) #如果使用ensg担心可能对不上\
  #symbols <- unique(Potential_DEG_8027$symbol)
  #outfile <- "DEGs\\test_promoter.bed"
  #type <- "bigWigAverageOverBed" #type -- "deeptools" has strand #"bed" only three col
  GTEx_gene_tab <- read.table("E:\\IDH_Hyper\\Data\\annotation\\GTEx_v19_gene_tab_ucsc.txt",
                              header=F, sep="\t", fill=T, check.names=F, quote="")
  
  sub_tab <- na.omit(GTEx_gene_tab[match(symbols, GTEx_gene_tab$V7), ])
  
  tss_bed <- c()
  for(i in 1:nrow(sub_tab)){
    tss <- 0
    if(sub_tab$V4[i] == "+"){
      tss <- as.numeric(sub_tab$V2[i])
    }
    else{
      tss <- as.numeric(sub_tab$V3[i])
    }
    promoter_start <- tss - upstream
    promoter_end <- tss + downstream
    if(promoter_start < 0){
      promoter_start <- 0
    }
    
    tss_bed <- rbind(tss_bed, c(promoter_start, promoter_end))
  }
  if(type == "bigWigAverageOverBed"){
    tss_bed4 <- data.frame(sub_tab$V1, tss_bed, paste(sub_tab$V7, sub_tab$V6, sep="<>") )
    sort_tss_bed4 <- tss_bed4[order(tss_bed4[,1], tss_bed4[,2]), ]
    write.table(sort_tss_bed4, outfile, sep="\t", row.names=F, col.names=F, quote=F)
  }
  if(type == "deeptools"){
    tss_bed6 <- data.frame(sub_tab$V1, tss_bed, sub_tab$V6, sub_tab$V7, sub_tab$V4 )
    sort_tss_bed6 <- tss_bed6[order(tss_bed6[,1], tss_bed6[,2]), ]
    write.table(sort_tss_bed6, outfile, sep="\t", row.names=F, col.names=F, quote=F)
  }
  if(type == "bed"){
    tss_bed3 <- data.frame(sub_tab$V1, tss_bed)
    sort_tss_bed3 <- tss_bed3[order(tss_bed3[,1], tss_bed3[,2]), ]
    write.table(sort_tss_bed3, outfile, sep="\t", row.names=F, col.names=F, quote=F)
  }
}
