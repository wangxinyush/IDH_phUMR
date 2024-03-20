#exp and methy - DEG (volcano)

setwd("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202203_exp_methy_plot")

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255) 
orange = rgb(243, 146, 0, max=255)

##Data1 - partially/fully hyper genes (symbol)
promoter_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\promoter_correct_refumr_anno_type.bed",
                             header=T, sep="\t", fill=T, check.names=F, quote="")
body_hyper <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\202201_promoter_correct\\result\\body_correct_refumr_anno_type.bed",
                         header=T, sep="\t", fill=T, check.names=F, quote="")

partial_hyper <- rbind(promoter_hyper[which(promoter_hyper$Type == "Patially"), ], body_hyper[which(body_hyper$Type == "Patially"), ])
full_hyper <- rbind(promoter_hyper[which(promoter_hyper$Type == "Fully"), ], body_hyper[which(body_hyper$Type == "Fully"), ])

partial_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Patially"), "Symbol"]), 
                          as.character(body_hyper[which(body_hyper$Type == "Patially"), "Symbol"])) )
full_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Fully"), "Symbol"]), 
                       as.character(body_hyper[which(body_hyper$Type == "Fully"), "Symbol"])) )

##Data2 - gene exp
TCGA_DEG <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202112_DEG\\TCGA_IDH_DEG_padjust.txt", 
                       header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")
DKFZ_DEG <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202112_DEG\\DKFZ_IDH_DEG_padjust.txt", 
                       header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

##Plot - volcano
###1. partially genes - DEG volcano
partial_plot_genes <-  unique(intersect(partial_hyper$Symbol, DKFZ_DEG$symbol))
partial_plot_genes_matrix <- c() #fc, padjust
for(i in 1:length(partial_plot_genes)){
  gene <- partial_plot_genes[i]
  partial_plot_genes_matrix <- rbind(partial_plot_genes_matrix, 
                             cbind( partial_hyper[which(partial_hyper$Symbol == gene), c(1:3, 6, 10:11, 17:19)], 
                             DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}
partial_plot_genes_matrix <- partial_plot_genes_matrix[which(partial_plot_genes_matrix[,"hyper_mean_IDH"] - partial_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
partial_up_index <- which(partial_plot_genes_matrix$fc_IDHn >= 2 & partial_plot_genes_matrix$padjust_IDHn <= 0.05)
partial_down_index <- which(partial_plot_genes_matrix$fc_IDHn <= 0.5 & partial_plot_genes_matrix$padjust_IDHn <= 0.05)

plot(log2(partial_plot_genes_matrix[,"fc_IDHn"]), -log10(partial_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(partial_plot_genes_matrix[partial_up_index,"fc_IDHn"]), -log10(partial_plot_genes_matrix[partial_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(partial_plot_genes_matrix[partial_down_index,"fc_IDHn"]), -log10(partial_plot_genes_matrix[partial_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###2. fully genes - DEG volcano
full_plot_genes <- unique(intersect( setdiff(full_hyper$Symbol, partial_hyper$Symbol), DKFZ_DEG$symbol))
full_plot_genes <- setdiff(full_plot_genes, c("KIF11", "IL17RD", "IRX2")) #three wrong genes
full_plot_genes_matrix <- c()
for(i in 1:length(full_plot_genes)){
  gene <- full_plot_genes[i]
  full_plot_genes_matrix <- rbind(full_plot_genes_matrix, 
                                  cbind( full_hyper[which(full_hyper$Symbol == gene), c(1:3, 6,10:11, 17:19)], 
                                         DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}
full_plot_genes_matrix <- full_plot_genes_matrix[which(full_plot_genes_matrix[,"hyper_mean_IDH"] - full_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
full_up_index <- which(full_plot_genes_matrix$fc_IDHn >= 2 & full_plot_genes_matrix$padjust_IDHn <= 0.05 )
full_down_index <- which(full_plot_genes_matrix$fc_IDHn <= 0.5 & full_plot_genes_matrix$padjust_IDHn <= 0.05 )

plot(log2(full_plot_genes_matrix[,"fc_IDHn"]), -log10(full_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(full_plot_genes_matrix[full_up_index,"fc_IDHn"]), -log10(full_plot_genes_matrix[full_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(full_plot_genes_matrix[full_down_index,"fc_IDHn"]), -log10(full_plot_genes_matrix[full_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###3. partially genes - DEG volcano (only promoter)
promoter_partial_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Patially"), "Symbol"])) )

promoter_partial_plot_genes <- unique(intersect(promoter_partial_genes, DKFZ_DEG$symbol))
promoter_partial_plot_genes_matrix <- c()
for(i in 1:length(promoter_partial_plot_genes)){
  gene <- promoter_partial_plot_genes[i]
  promoter_partial_plot_genes_matrix <- rbind(promoter_partial_plot_genes_matrix, 
                                              cbind( partial_hyper[which(partial_hyper$Symbol == gene), c(1:3, 6,10:11, 17:19)], 
                                                     DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}
promoter_partial_plot_genes_matrix <- promoter_partial_plot_genes_matrix[which(promoter_partial_plot_genes_matrix[,"hyper_mean_IDH"] - promoter_partial_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
promoter_partial_up_index <- which(promoter_partial_plot_genes_matrix$fc_IDHn >= 2 & promoter_partial_plot_genes_matrix$padjust_IDHn <= 0.05 )
promoter_partial_down_index <- which(promoter_partial_plot_genes_matrix$fc_IDHn <= 0.5 & promoter_partial_plot_genes_matrix$padjust_IDHn <= 0.05 )

plot(log2(promoter_partial_plot_genes_matrix[,"fc_IDHn"]), -log10(promoter_partial_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(promoter_partial_plot_genes_matrix[promoter_partial_up_index,"fc_IDHn"]), -log10(promoter_partial_plot_genes_matrix[promoter_partial_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(promoter_partial_plot_genes_matrix[promoter_partial_down_index,"fc_IDHn"]), -log10(promoter_partial_plot_genes_matrix[promoter_partial_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###4. fully genes - DEG volcano (only promoter)
promoter_full_genes <- unique(c(as.character(promoter_hyper[which(promoter_hyper$Type == "Fully"), "Symbol"])) )
promoter_full_plot_genes <- unique(intersect( setdiff(promoter_full_genes, promoter_partial_genes), DKFZ_DEG$symbol))
promoter_full_plot_genes <- setdiff(promoter_full_plot_genes, c("KIF11", "IL17RD", "IRX2")) #three wrong genes
promoter_full_plot_genes_matrix <- c()
for(i in 1:length(promoter_full_plot_genes)){
  gene <- promoter_full_plot_genes[i]
  promoter_full_plot_genes_matrix <- rbind(promoter_full_plot_genes_matrix, 
                                           cbind( promoter_hyper[which(promoter_hyper$Symbol == gene), c(1:3, 6,10:11, 17:19)], 
                                                  DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}

promoter_full_plot_genes_matrix <- promoter_full_plot_genes_matrix[which(promoter_full_plot_genes_matrix[,"hyper_mean_IDH"] - promoter_full_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
promoter_full_up_index <- which(promoter_full_plot_genes_matrix$fc_IDHn >= 2 & promoter_full_plot_genes_matrix$padjust_IDHn <= 0.05 )
#promoter_full_plot_genes_matrix[promoter_full_up_index,] #check 8 genes
promoter_full_down_index <- which(promoter_full_plot_genes_matrix$fc_IDHn <= 0.5 & promoter_full_plot_genes_matrix$padjust_IDHn <= 0.05 )

plot(log2(promoter_full_plot_genes_matrix[,"fc_IDHn"]), -log10(promoter_full_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(promoter_full_plot_genes_matrix[promoter_full_up_index,"fc_IDHn"]), -log10(promoter_full_plot_genes_matrix[promoter_full_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(promoter_full_plot_genes_matrix[promoter_full_down_index,"fc_IDHn"]), -log10(promoter_full_plot_genes_matrix[promoter_full_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###5. partially genes - DEG volcano (only promoter)
body_partial_genes <- unique(c(as.character(body_hyper[which(body_hyper$Type == "Patially"), "Symbol"])) )

body_partial_plot_genes <- unique(intersect(body_partial_genes, DKFZ_DEG$symbol))
body_partial_plot_genes_matrix <- c()
for(i in 1:length(body_partial_plot_genes)){
  gene <- body_partial_plot_genes[i]
  body_partial_plot_genes_matrix <- rbind(body_partial_plot_genes_matrix, 
                                          cbind( partial_hyper[which(partial_hyper$Symbol == gene), c(1:3, 6,10:11, 17:19)], 
                                                 DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}
body_partial_plot_genes_matrix <- body_partial_plot_genes_matrix[which(body_partial_plot_genes_matrix[,"hyper_mean_IDH"] - body_partial_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
body_partial_up_index <- which(body_partial_plot_genes_matrix$fc_IDHn >= 2 & body_partial_plot_genes_matrix$padjust_IDHn <= 0.05 )
body_partial_down_index <- which(body_partial_plot_genes_matrix$fc_IDHn <= 0.5 & body_partial_plot_genes_matrix$padjust_IDHn <= 0.05 )

plot(log2(body_partial_plot_genes_matrix[,"fc_IDHn"]), -log10(body_partial_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(body_partial_plot_genes_matrix[body_partial_up_index,"fc_IDHn"]), -log10(body_partial_plot_genes_matrix[body_partial_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(body_partial_plot_genes_matrix[body_partial_down_index,"fc_IDHn"]), -log10(body_partial_plot_genes_matrix[body_partial_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###6. fully genes - DEG volcano (only promoter)
body_full_genes <- unique(c(as.character(body_hyper[which(body_hyper$Type == "Fully"), "Symbol"])) )
body_full_plot_genes <- unique(intersect( setdiff(body_full_genes, body_partial_genes), DKFZ_DEG$symbol))
body_full_plot_genes_matrix <- c()
for(i in 1:length(body_full_plot_genes)){
  gene <- body_full_plot_genes[i]
  body_full_plot_genes_matrix <- rbind(body_full_plot_genes_matrix, 
                                       cbind( body_hyper[which(body_hyper$Symbol == gene), c(1:3, 6,10:11, 17:19)], 
                                              DKFZ_DEG[which(DKFZ_DEG$symbol == gene), c(1:6, 9,10)] ) )
}

body_full_plot_genes_matrix <- body_full_plot_genes_matrix[which(body_full_plot_genes_matrix[,"hyper_mean_IDH"] - body_full_plot_genes_matrix[,"hyper_mean_normal"] >= 0.2), ]
body_full_up_index <- which(body_full_plot_genes_matrix$fc_IDHn >= 2 & body_full_plot_genes_matrix$padjust_IDHn <= 0.05 )
#body_full_plot_genes_matrix[body_full_up_index,] #check 8 genes
body_full_down_index <- which(body_full_plot_genes_matrix$fc_IDHn <= 0.5 & body_full_plot_genes_matrix$padjust_IDHn <= 0.05 )

plot(log2(body_full_plot_genes_matrix[,"fc_IDHn"]), -log10(body_full_plot_genes_matrix[,"padjust_IDHn"]), 
     main = "IDH vs. Normal(DKFZ)", xlab="log2(fold change)", ylab="-log10(p.adjust)", 
     pch=20, col="gray", cex.lab=1.3, las=1, bty="l")

points(log2(body_full_plot_genes_matrix[body_full_up_index,"fc_IDHn"]), -log10(body_full_plot_genes_matrix[body_full_up_index,"padjust_IDHn"]),
       pch=20, col=red)
points(log2(body_full_plot_genes_matrix[body_full_down_index,"fc_IDHn"]), -log10(body_full_plot_genes_matrix[body_full_down_index,"padjust_IDHn"]),
       pch=20, col=blue)

###4. Check down DEG - 9 genes
promoter_full_plot_genes_matrix[promoter_full_up_index,]



###write.table
write.table( partial_plot_genes_matrix, "hyper_gene_matrix\\partial_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( partial_plot_genes_matrix[partial_up_index, ], "hyper_gene_matrix\\partial_up_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( partial_plot_genes_matrix[partial_down_index, ], "hyper_gene_matrix\\partial_down_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )

write.table( full_plot_genes_matrix, "hyper_gene_matrix\\full_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( full_plot_genes_matrix[full_up_index, ], "hyper_gene_matrix\\full_up_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( full_plot_genes_matrix[full_down_index, ], "hyper_gene_matrix\\full_down_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )

write.table( promoter_partial_plot_genes_matrix, "hyper_gene_matrix\\promoter_partial_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( promoter_partial_plot_genes_matrix[promoter_partial_up_index, ], "hyper_gene_matrix\\promoter_partial_up_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( promoter_partial_plot_genes_matrix[promoter_partial_down_index, ], "hyper_gene_matrix\\promoter_partial_down_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )

write.table( promoter_full_plot_genes_matrix, "hyper_gene_matrix\\promoter_full_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( promoter_full_plot_genes_matrix[promoter_full_up_index, ], "hyper_gene_matrix\\promoter_full_up_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( promoter_full_plot_genes_matrix[promoter_full_down_index, ], "hyper_gene_matrix\\promoter_full_down_genes_matrix.txt", sep="\t", row.names=F, col.names=T, quote=F )


###write.table - partial up/down hyper regions
