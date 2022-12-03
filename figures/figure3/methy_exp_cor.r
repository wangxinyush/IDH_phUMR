#@2022.05.11
#correlation of methylation changes and expression

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255) 
orange = rgb(243, 146, 0, max=255)

#1. Data preparasion
#Data1 - partially/fully hyper genes (symbol)
promoter_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\promoter_correct_refumr_anno_type.bed",
                             header=T, sep="\t", fill=T, check.names=F, quote="")
body_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\body_correct_refumr_anno_type.bed",
                         header=T, sep="\t", fill=T, check.names=F, quote="")

hyper_type_position_anno <- rbind( cbind( promoter_hyper, postion=rep("Promoter", nrow(promoter_hyper))), 
       cbind( body_hyper, postion=rep("Body", nrow(body_hyper))) ) #partially/fully, promoter/body

#because of repetitive
#which(hyper_type_position_anno$hyper_start == 7449725) #147 1769
#which(hyper_type_position_anno$hyper_start == 77159475) #182 1784
#which(hyper_type_position_anno$hyper_start == 71630574) #1134 2093

hyper_type_position_anno <- hyper_type_position_anno[-c(1769, 1784, 2093),]
rownames(hyper_type_position_anno) <- paste(hyper_type_position_anno$chr, hyper_type_position_anno$hyper_start, hyper_type_position_anno$hyper_end, sep=",")

#Data2 - methylation
partial_methy <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202204_annotation_three_regions\\methy_matrix\\sort_partial_hyper_methy_matrix.txt",
                            header=T, sep="\t", fill=T, check.names=F, quote="")
full_methy <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202204_annotation_three_regions\\methy_matrix\\sort_full_hyper_methy_matrix.txt",
                         header=T, sep="\t", fill=T, check.names=F, quote="")
hyper_methy <- rbind(partial_methy, full_methy)
rownames(hyper_methy) <- paste(hyper_methy$chr, hyper_methy$start, hyper_methy$end, sep=",")

DKFZ_hyper_methy <- hyper_methy[, 79:138]
colname_arr <- c()
for(i in 1:ncol(DKFZ_hyper_methy)){
  colname_arr <- c(colname_arr, strsplit(colnames(DKFZ_hyper_methy)[i], "_")[[1]][2])
}
colnames(DKFZ_hyper_methy) <- colname_arr

#Data3 - exp
DKFZ_gbm_tpm <- read.table("E:\\IDH_Hyper\\Data\\60_GBM_Data\\GSE121720_RNAseq_expression_matrix_TPMs.txt", header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

col_index <- match( colnames(DKFZ_hyper_methy), colnames(DKFZ_gbm_tpm) )
colnames(DKFZ_gbm_tpm)[col_index]  #check right

#anno
DKFZ_GEN_anno <- read.table("E:\\IDH_Hyper\\Data\\60_GBM_Data\\gencode.v19.gene.map",
                            header=F, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

#2. Correlation matrix
#region, symbol, ensg, cor_pearson, cor_spearman
methy_exp_cor <- c()
for (i in 1:nrow(hyper_type_position_anno)){ #i = 76, SYT6
  region <- rownames(hyper_type_position_anno)[i]
  symbol <- as.character(hyper_type_position_anno$Symbol[i])
  ensg <- rownames(DKFZ_GEN_anno)[ which(DKFZ_GEN_anno$V2 == symbol)[1] ]
  
  type <- as.character(hyper_type_position_anno$Type[i])
  position <- as.character(hyper_type_position_anno$postion[i])
  
  methy <- as.numeric(DKFZ_hyper_methy[region,])
  exp <- as.numeric(DKFZ_gbm_tpm[ensg, col_index])
  
  cor_pearson <- round( cor(methy, exp, method="pearson"), 3)
  cor_spearman <- round( cor(methy, exp, method="spearman"), 3)
  
  methy_exp_cor <- rbind(methy_exp_cor, c(region, symbol, ensg, type, position, cor_pearson, cor_spearman))
}

colnames(methy_exp_cor) <- c("region", "symbol", "ensg", "type", "position", "cor_pearson", "cor_spearman")

partial_index <- which(methy_exp_cor[,"type"] == "Patially")
full_index <- which(methy_exp_cor[,"type"] == "Fully")

all_promoter_index <- which(methy_exp_cor[,"position"] == "Promoter")
partial_promoter_index <- which(methy_exp_cor[,"type"] == "Patially" & methy_exp_cor[,"position"] == "Promoter")
full_promoter_index <- which(methy_exp_cor[,"type"] == "Fully" & methy_exp_cor[,"position"] == "Promoter")

all_body_index <- which(methy_exp_cor[,"position"] == "Body")
partial_body_index <- which(methy_exp_cor[,"type"] == "Patially" & methy_exp_cor[,"position"] == "Body")
full_body_index <- which(methy_exp_cor[,"type"] == "Fully" & methy_exp_cor[,"position"] == "Body")

#all - 9
boxplot(na.omit(as.numeric(methy_exp_cor[, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[partial_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[full_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_pearson"])),
        col=c("gray", red, orange), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 1, outline = F, ylim = c(-1, 1),
        #names=c("All", "Partially", "Fully"),
        main="All/Promoter/Body correlation" )
dev.print(pdf, file="exp_methy_Rplot/all_promoter_cor.pdf", width=6, height=8)

#all
t.test( na.omit(as.numeric(methy_exp_cor[, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_index, "cor_pearson"])) )$p.value #0.670
t.test( na.omit(as.numeric(methy_exp_cor[, "cor_pearson"])) ,  
       na.omit(as.numeric(methy_exp_cor[full_index, "cor_pearson"])) )$p.value #0.504191
t.test( na.omit(as.numeric(methy_exp_cor[partial_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_index, "cor_pearson"])) )$p.value #0.3661653

#promoter
t.test( na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_pearson"])) )$p.value #0.01145379
t.test( na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_pearson"])) )$p.value #7.895532e-08
t.test( na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_pearson"])) )$p.value #7.8265e-11

#body
t.test( na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_pearson"])) )$p.value #0.04082996
t.test( na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_pearson"])) )$p.value #0.09196358
t.test( na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_pearson"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_pearson"])) )$p.value #0.001257089

#promoter
boxplot(na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[parital_promoter_index, "cor_pearson"])),
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_pearson"])),
        col=c("gray", red, orange), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 1, outline = F, ylim = c(-1, 1),
        names=c("All", "Partially", "Fully"),
        main="Promoter correlation" )
dev.print(pdf, file="exp_methy_Rplot/promoter_cor.pdf", width=6, height=8)



###Spearman
boxplot(na.omit(as.numeric(methy_exp_cor[, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[partial_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[full_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_spearman"])),
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_spearman"])),
        col=c("gray", red, orange), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 1, outline = F, ylim = c(-1, 1),
        #names=c("All", "Partially", "Fully"),
        main="All/Promoter/Body correlation" )
dev.print(pdf, file="exp_methy_Rplot/all_promoter_cor_spearman.pdf", width=6, height=8)

#all
t.test( na.omit(as.numeric(methy_exp_cor[, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_index, "cor_spearman"])) )$p.value #0.4898843
t.test( na.omit(as.numeric(methy_exp_cor[, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_index, "cor_spearman"])) )$p.value #0.2751577
t.test( na.omit(as.numeric(methy_exp_cor[partial_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_index, "cor_spearman"])) )$p.value #0.1402333

#promoter
t.test( na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_spearman"])) )$p.value #0.004763874
t.test( na.omit(as.numeric(methy_exp_cor[all_promoter_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_spearman"])) )$p.value #2.423088e-09
t.test( na.omit(as.numeric(methy_exp_cor[partial_promoter_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_promoter_index, "cor_spearman"])) )$p.value #5.336301e-13

#body
t.test( na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_spearman"])) )$p.value #0.07551778
t.test( na.omit(as.numeric(methy_exp_cor[all_body_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_spearman"])) )$p.value #0.123131
t.test( na.omit(as.numeric(methy_exp_cor[partial_body_index, "cor_spearman"])) ,  
        na.omit(as.numeric(methy_exp_cor[full_body_index, "cor_spearman"])) )$p.value #0.004041299
