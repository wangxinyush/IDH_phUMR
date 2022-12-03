#@2021.12.26
#DEGs of GTEx/TCGA

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255)
orange = rgb(243, 146, 0, max=255)

####1. Data & Samples
#TPM Matrix - GTEx & TCGA
Brain_GTEx_TCGA_tpm <- read.table(gzfile("E:\\IDH_Hyper\\Data\\GTEx_TCGA_Exp\\Brain_GTEx_TCGA_gene_tpm.gz"), 
                                  header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")
Brain_GTEx_TCGA_tpm <- log2(2^Brain_GTEx_TCGA_tpm - 0.001 + 1) #raw log2(tpm+0.001) is not suitable for visualization 

#samples & gene annotation
#TCGA_Exp_samples - 697 samples
#IDH:434
TCGA_exp_IDH_samples <- read.table("E:\\IDH_Hyper\\Results\\3_Oncogenic_pathway\\TCGA_samples\\TCGA_exp_IDH_samples.txt",
                                   header=F, sep="\t", fill=T, check.names=F, quote="")
#WT:263
TCGA_exp_WT_samples <- read.table("E:\\IDH_Hyper\\Results\\3_Oncogenic_pathway\\TCGA_samples\\TCGA_exp_WT_samples.txt",
                                  header=F, sep="\t", fill=T, check.names=F, quote="")



#GTEx samples - 2076
GTEx_Brain_samples <- read.table("E:\\IDH_Hyper\\Data\\GTEx_TCGA_Exp\\GTEx_Brain_samples.txt",
                                 header=F, sep="\t", fill=T, check.names=F, quote="")
#Gene Annotation
#60498 x 5
GEN_anno <- read.table("E:\\IDH_Hyper\\Data\\GTEx_TCGA_Exp\\gencode.v23.annotation.gene.probemap",
                       header=T, row.names = 1, sep="\t", fill=T, check.names=F, quote="")

#Aviable Normal, IDH, WT samples
#Brain: 1146, IDH: 427, WT: 261
IDH_samples <- intersect(colnames(Brain_GTEx_TCGA_tpm), as.character(TCGA_exp_IDH_samples$V1))
WT_samples <- intersect(colnames(Brain_GTEx_TCGA_tpm), as.character(TCGA_exp_WT_samples$V1))
Brain_samples <- intersect(colnames(Brain_GTEx_TCGA_tpm), as.character(GTEx_Brain_samples$V1))


##Check pheatmap for partially up-regulated genes @2022.11.3
partial_DKFZ_up_genes <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\partial_DEG_genes\\partial_DKFZ_up_genes.txt",
                       header=F, sep="\t", fill=T, check.names=F, quote="")
partial_DKFZ_down_genes <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\partial_DEG_genes\\partial_DKFZ_down_genes.txt",
                                    header=F, sep="\t", fill=T, check.names=F, quote="")
partial_TCGA_up_genes <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\partial_DEG_genes\\partial_TCGA_up_genes.txt",
                                    header=F, sep="\t", fill=T, check.names=F, quote="")
partial_TCGA_down_genes <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\partial_DEG_genes\\partial_TCGA_down_genes.txt",
                                    header=F, sep="\t", fill=T, check.names=F, quote="")

library(pheatmap)
pheatmap(Brain_GTEx_TCGA_tpm[match( partial_DKFZ_up_genes$V1, GEN_anno$gene), Brain_samples],
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = F, show_colnames = F, cluster_cols = T,
         clustering_method = "ward.D2")
pheatmap(Brain_GTEx_TCGA_tpm[match( partial_DKFZ_down_genes$V1, GEN_anno$gene), Brain_samples],
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = F, show_colnames = F, cluster_cols = T,
         clustering_method = "ward.D2")
pheatmap(Brain_GTEx_TCGA_tpm[match( partial_TCGA_up_genes$V1, GEN_anno$gene), Brain_samples],
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = F, show_colnames = F, cluster_cols = T,
         clustering_method = "ward.D2")
pheatmap(Brain_GTEx_TCGA_tpm[match( partial_TCGA_down_genes$V1, GEN_anno$gene), Brain_samples],
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = F, show_colnames = F, cluster_cols = T,
         clustering_method = "ward.D2")


####2. DEG (IDH vs. Normal, IDH vs. WT) - GTEx/TCGA
IDH_Normal_WT_DEG <- c()
IDH_Normal_WT_DEG_rownames <- c()
IDH_Normal_WT_DEG_genes <- c() #ensg(name): symbol(value)
#for(i in 1:100){
for(i in 1:nrow(Brain_GTEx_TCGA_tpm)){
  ensg <- as.character(rownames(Brain_GTEx_TCGA_tpm)[i])
  Normal_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, Brain_samples]))
  IDH_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, IDH_samples]))
  WT_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, WT_samples]))
  
  mean_Normal <- mean(Normal_exp)
  mean_IDH <- mean(IDH_exp)
  mean_WT <- mean(WT_exp)
  
  fc_IDHn <- 2^(mean_IDH - mean_Normal) #tpm is log2(tpm + 1) #IDH vs. Normal
  fc_IDHw <- 2^(mean_IDH - mean_WT) #IDH vs. WT
  
  pvalue_IDHn <- 1
  pvalue_IDHw <- 1
  #pvalue: IDH vs. Normal
  if(length(IDH_exp) >= 3 & length(Normal_exp) >= 3 & mean_IDH - mean_Normal != 0){
    pvalue_IDHn <- t.test(IDH_exp, Normal_exp)$p.value
  }
  #pvalue: IDH vs. WT
  if(length(IDH_exp) >= 3 & length(WT_exp) >= 3 & mean_IDH - mean_WT != 0){
    pvalue_IDHw <- t.test(IDH_exp, WT_exp)$p.value
  }
  
  IDH_Normal_WT_DEG <- rbind(IDH_Normal_WT_DEG, c(mean_Normal, mean_IDH, mean_WT, fc_IDHn, fc_IDHw, pvalue_IDHn, pvalue_IDHw))
  IDH_Normal_WT_DEG_rownames <- c(IDH_Normal_WT_DEG_rownames, ensg)
  IDH_Normal_WT_DEG_genes <- c(IDH_Normal_WT_DEG_genes, as.character(GEN_anno[ensg,"gene"]))
}
rownames(IDH_Normal_WT_DEG) <- IDH_Normal_WT_DEG_rownames
names(IDH_Normal_WT_DEG_genes) <- IDH_Normal_WT_DEG_rownames
padjust_IDHn <- p.adjust(IDH_Normal_WT_DEG[,6], method = "BH", n = nrow(IDH_Normal_WT_DEG))
padjust_IDHw <- p.adjust(IDH_Normal_WT_DEG[,7], method = "BH", n = nrow(IDH_Normal_WT_DEG))

IDH_Normal_WT_DEG_padjust <- data.frame(IDH_Normal_WT_DEG_genes, IDH_Normal_WT_DEG, padjust_IDHn, padjust_IDHw)
colnames(IDH_Normal_WT_DEG_padjust) <- c("symbol", "mean_Normal_log2", "mean_IDH_log2", "mean_WT_log2", "fc_IDHn", "fc_IDHw",
                              "pvalue_IDHn", "pvalue_IDHw", "padjust_IDHn", "padjust_IDHw")

write.table( IDH_Normal_WT_DEG_padjust, "TCGA_IDH_DEG_padjust.txt", sep="\t", row.names=T, col.names=T, quote=F )


##check: ENSG00000134207.14	SYT6, ENSG00000159713.10 TPPP3
ensg <- "ENSG00000134207.14" #SYT6 - right
ensg <- "ENSG00000159713.10" #TPPP3 - right
Normal_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, Brain_samples]))
IDH_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, IDH_samples]))
WT_exp <- na.omit(as.numeric(Brain_GTEx_TCGA_tpm[ensg, WT_samples]))

mean(Normal_exp)
mean(IDH_exp)
mean(WT_exp)
2^(mean(IDH_exp) - mean(Normal_exp))
2^(mean(IDH_exp) - mean(WT_exp))
t.test(IDH_exp, Normal_exp)$p.value
t.test(IDH_exp, WT_exp)$p.value
