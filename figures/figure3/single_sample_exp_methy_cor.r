#correlation of expression (log2FC) and methy changes in AK213 and AK071 

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255) 
orange = rgb(243, 146, 0, max=255)


##1) Data1 - partially/fully hyper genes (symbol)
promoter_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\promoter_correct_refumr_anno_type.bed",
                             header=T, sep="\t", fill=T, check.names=F, quote="")
body_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\body_correct_refumr_anno_type.bed",
                         header=T, sep="\t", fill=T, check.names=F, quote="")

partial_hyper <- rbind(promoter_hyper[which(promoter_hyper$Type == "Patially"), ], body_hyper[which(body_hyper$Type == "Patially"), ])
full_hyper <- rbind(promoter_hyper[which(promoter_hyper$Type == "Fully"), ], body_hyper[which(body_hyper$Type == "Fully"), ])

#rownames(partial_hyper) <- paste(partial_hyper$chr, partial_hyper$hyper_start, partial_hyper$hyper_end, sep=",")
#rownames(full_hyper) <- paste(full_hyper$chr, full_hyper$hyper_start, full_hyper$hyper_end, sep=",")

partial_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Patially"), "Symbol"]), 
                          as.character(body_hyper[which(body_hyper$Type == "Patially"), "Symbol"])) )
full_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Fully"), "Symbol"]), 
                       as.character(body_hyper[which(body_hyper$Type == "Fully"), "Symbol"])) )

#parital methy
partial_methy <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202204_annotation_three_regions\\methy_matrix\\sort_partial_hyper_methy_matrix.txt",
                            header=T, sep="\t", fill=T, check.names=F, quote="")
rownames(partial_methy) <- paste(partial_methy$chr, partial_methy$start, partial_methy$end, sep=",")

partial_methy_normal_IDH_WT <- c() #normal, AK213, AK071
for(i in 1:nrow(partial_methy)){
  normal_methy <- na.omit(as.numeric(partial_methy[i,c(4:6)]))
  partial_methy_normal_IDH_WT <- rbind(partial_methy_normal_IDH_WT, 
                                     c(mean(normal_methy), partial_methy[i, "GSM3444668_AK213"], partial_methy[i, "GSM3444631_AK071"]) )
}
rownames(partial_methy_normal_IDH_WT) <- rownames(partial_methy)
colnames(partial_methy_normal_IDH_WT) <- c("normal", "AK213", "AK071")

#full_methy
full_methy <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202204_annotation_three_regions\\methy_matrix\\sort_full_hyper_methy_matrix.txt",
                            header=T, sep="\t", fill=T, check.names=F, quote="")
rownames(full_methy) <- paste(full_methy$chr, full_methy$start, full_methy$end, sep=",")

full_methy_normal_IDH_WT <- c() #normal, AK213, AK071
for(i in 1:nrow(full_methy)){
  normal_methy <- na.omit(as.numeric(full_methy[i,c(4:6)]))
  full_methy_normal_IDH_WT <- rbind(full_methy_normal_IDH_WT, 
                                       c(mean(normal_methy), full_methy[i, "GSM3444668_AK213"], full_methy[i, "GSM3444631_AK071"]) )
}
rownames(full_methy_normal_IDH_WT) <- rownames(full_methy)
colnames(full_methy_normal_IDH_WT) <- c("normal", "AK213", "AK071")


##2) Data2 - Gene Expression
DKFZ_gbm_tpm <- read.table("E:\\IDH_Hyper\\Data\\60_GBM_Data\\GSE121720_RNAseq_expression_matrix_TPMs.txt", header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

#normal mean
DKFZ_normal_index <- c(71:75)
DKFZ_exp_normal_IDH_WT <- c()
for(i in 1:nrow(DKFZ_gbm_tpm)){
  normal_exp <- na.omit(as.numeric(DKFZ_gbm_tpm[i,DKFZ_normal_index]))
  DKFZ_exp_normal_IDH_WT <- rbind(DKFZ_exp_normal_IDH_WT, 
                           c(mean(normal_exp), DKFZ_gbm_tpm[i, "AK213"], DKFZ_gbm_tpm[i, "AK071"]) )
}
rownames(DKFZ_exp_normal_IDH_WT) <- rownames(DKFZ_gbm_tpm)
colnames(DKFZ_exp_normal_IDH_WT) <- c("normal", "AK213", "AK071")

#anno
DKFZ_GEN_anno <- read.table("E:\\IDH_Hyper\\Data\\60_GBM_Data\\gencode.v19.gene.map",
                            header=F, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

##3) methy_exp matrix
#which(partial_hyper$Symbol == "MET")[1] #three regions
partial_gene_methy_exp <- c() #symbol, ensg, region, methy_AK213, methy_normal, exp_AK213, exp_normal.
for(i in 1:length(partial_genes)){
  symbol <- partial_genes[i]
  region_i <- which(partial_hyper$Symbol == symbol)[1]
  
  hyper_region <- paste(partial_hyper$chr[region_i], partial_hyper$hyper_start[region_i], partial_hyper$hyper_end[region_i], sep=",")
  ensg <- rownames(DKFZ_GEN_anno)[ which(DKFZ_GEN_anno$V2 == symbol)[1] ]
  
  if(!is.na(ensg)){
    partial_gene_methy_exp <- rbind(partial_gene_methy_exp, 
                                    c(symbol, ensg, hyper_region,
                                      partial_methy_normal_IDH_WT[hyper_region, ],
                                      DKFZ_exp_normal_IDH_WT[ensg, ]
                                    ) )
  }
}
colnames(partial_gene_methy_exp) <- c("symbol", "ensg", "region", "methy_normal", "methy_IDH", "methy_WT",
                                      "exp_normal", "exp_IDH", "exp_WT")
full_gene_methy_exp <- c() #symbol, ensg, region, methy_AK213, methy_normal, exp_AK213, exp_normal.
for(i in 1:length(full_genes)){
  symbol <- full_genes[i]
  region_i <- which(full_hyper$Symbol == symbol)[1]
  
  hyper_region <- paste(full_hyper$chr[region_i], full_hyper$hyper_start[region_i], full_hyper$hyper_end[region_i], sep=",")
  ensg <- rownames(DKFZ_GEN_anno)[ which(DKFZ_GEN_anno$V2 == symbol)[1] ]
  
  if(!is.na(ensg)){
    full_gene_methy_exp <- rbind(full_gene_methy_exp, 
                                 c(symbol, ensg, hyper_region,
                                   full_methy_normal_IDH_WT[hyper_region, ],
                                   DKFZ_exp_normal_IDH_WT[ensg, ]
                                 ) )
  }
}
colnames(full_gene_methy_exp) <- c("symbol", "ensg", "region", "methy_normal", "methy_IDH", "methy_WT",
                                      "exp_normal", "exp_IDH", "exp_WT")

#check right: - HES5
partial_hyper[2,] #HES5, chr1 2460471   2461446
head(partial_methy_normal_IDH_WT, 2) #0.02333333 0.789 0.099
rownames(DKFZ_GEN_anno)[173] #ENSG00000197921.5	HES5
DKFZ_exp_normal_IDH_WT["ENSG00000197921.5",] #1.4670030 16.5604325  0.8475512
partial_gene_methy_exp[2,]

##4) Correlation (for a single sample)
#1. Data Preparation
#partial_gene_methy_exp - remove NA and to 1
normal_na_index <- which(is.na(partial_gene_methy_exp[,"methy_normal"]))
IDH_na_index <- which(is.na(partial_gene_methy_exp[,"methy_IDH"]))
WT_na_index <- which(is.na(partial_gene_methy_exp[,"methy_WT"]))

after_rmNA_index <- setdiff(1:nrow(partial_gene_methy_exp), unique(c(normal_na_index, IDH_na_index, WT_na_index)))

partial_gene_methy_exp_rmNA <- partial_gene_methy_exp[after_rmNA_index,]

#full_gene_methy_exp - remove NA and to 1
normal_na_index <- which(is.na(full_gene_methy_exp[,"methy_normal"]))
IDH_na_index <- which(is.na(full_gene_methy_exp[,"methy_IDH"]))
WT_na_index <- which(is.na(full_gene_methy_exp[,"methy_WT"]))

after_rmNA_index <- setdiff(1:nrow(full_gene_methy_exp), unique(c(normal_na_index, IDH_na_index, WT_na_index)))

full_gene_methy_exp_rmNA <- full_gene_methy_exp[after_rmNA_index,]

#2. correlation for a single sample
methy <- 


##4) bak - Correlation (delat methy, log2FC)
dim(partial_gene_methy_exp) #1481
dim(full_gene_methy_exp) #532

#1. partial - all , promoter, body
#partial_gene_methy_exp - remove NA and to 1
normal_na_index <- which(is.na(partial_gene_methy_exp[,"methy_normal"]))
IDH_na_index <- which(is.na(partial_gene_methy_exp[,"methy_IDH"]))
WT_na_index <- which(is.na(partial_gene_methy_exp[,"methy_WT"]))

after_rmNA_index <- setdiff(1:nrow(partial_gene_methy_exp), unique(c(normal_na_index, IDH_na_index, WT_na_index)))

partial_gene_methy_exp_rmNA <- partial_gene_methy_exp[after_rmNA_index,]

#all
all_partial_plot <- partial_gene_methy_exp_rmNA
change_methy <- as.numeric(all_partial_plot[,"methy_IDH"]) - as.numeric(all_partial_plot[,"methy_normal"])
fc_exp <- log2( (as.numeric(all_partial_plot[,"exp_IDH"]) + 1) / (as.numeric(all_partial_plot[,"exp_normal"]) + 1) )

correlation_methy_exp_plot(change_methy, fc_exp, "exp_methy_Rplot/all_exp_methy_cor", "IDH vs. Normal(DKFZ)" )

#promoter
all_full_plot <- full_gene_methy_exp_rmNA
single_methy <- as.numeric(all_full_plot[,"methy_IDH"])
log_exp <- log2( (as.numeric(all_full_plot[,"exp_IDH"]) + 0.1) )

correlation_methy_exp_plot(single_methy, log_exp, "exp_methy_Rplot/all_exp_methy_cor", "IDH vs. Normal(DKFZ)" )


#2. full - all , promoter, body
#all
#full_gene_methy_exp - remove NA and to 1
normal_na_index <- which(is.na(full_gene_methy_exp[,"methy_normal"]))
IDH_na_index <- which(is.na(full_gene_methy_exp[,"methy_IDH"]))
WT_na_index <- which(is.na(full_gene_methy_exp[,"methy_WT"]))

after_rmNA_index <- setdiff(1:nrow(full_gene_methy_exp), unique(c(normal_na_index, IDH_na_index, WT_na_index)))

full_gene_methy_exp_rmNA <- full_gene_methy_exp[after_rmNA_index,]

#all
all_full_plot <- full_gene_methy_exp_rmNA
change_methy <- as.numeric(all_full_plot[,"methy_IDH"]) - as.numeric(all_full_plot[,"methy_normal"])
fc_exp <- log2( (as.numeric(all_full_plot[,"exp_IDH"]) + 1) / (as.numeric(all_full_plot[,"exp_normal"]) + 1) )

correlation_methy_exp_plot(change_methy, fc_exp, "exp_methy_Rplot/all_exp_methy_cor", "IDH vs. Normal(DKFZ)" )

#promoter
promoter_full_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Fully"), "Symbol"]) ))
promoter_full_plot <- all_full_plot [ na.omit(match(promoter_full_genes, all_full_plot[,"symbol"])), ]

#length(intersect(promoter_full_plot[,"symbol"], promoter_full_genes))
change_methy <- as.numeric(promoter_full_plot[,"methy_IDH"]) - as.numeric(promoter_full_plot[,"methy_normal"])
fc_exp <- log2( (as.numeric(promoter_full_plot[,"exp_IDH"]) + 1) / (as.numeric(promoter_full_plot[,"exp_normal"]) + 1) )

correlation_methy_exp_plot(change_methy, fc_exp, "exp_methy_Rplot/promoter_full_exp_methy_cor", "IDH vs. Normal(DKFZ)" )



correlation_methy_exp_plot	<- function(change_methy, fc_exp, filename, main_name){
  cor_str <- paste(round( cor(change_methy, fc_exp, method="pearson"), 3),
                   round( cor(change_methy, fc_exp, method="spearman"), 3), sep=",")
  plot( change_methy, fc_exp,
        pch=20, col="gray", cex.lab=1.3, las=1, bty="l", font.axis=2, cex=2,
        main = "IDH vs. Normal(DKFZ) ", xlab=cor_str, ylab="log2(FC)",
        xlim=c(-0.2,1), ylim=c(-6,6) )
  abline(lm(fc_exp ~ change_methy), col = red, lwd=3)
  dev.print(png, file=paste(filename, ".png", sep=""), width=300, height=300)
  dev.print(pdf, file=paste(filename, ".pdf", sep=""), width=6, height=6)
  return(0)
}
