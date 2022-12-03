DEG <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202112_DEG\\DKFZ_IDH_DEG_padjust.txt",
                                  header=T, sep="\t", fill=T, check.names=F, quote="")

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255)
orange = rgb(243, 146, 0, max=255)
cgi_blue = rgb(57, 127, 185, max=255) #CGI

#取Cancer Pathway
pathway_genes <- read.table("E:\\IDH_Hyper\\Results\\3_Oncogenic_pathway\\pathway_mapper\\Cancer_Pathway_TCGA.txt",
                            header=T, sep="\t", fill=T, check.names=F, quote="")
length( as.character(pathway_genes[which(pathway_genes$`OG/TSG` == "OG"), 1]) ) #127

og_data <- read.table("E:\\IDH_Hyper\\Data\\COSMIC\\COSMIC_oncogene.txt",
                      header=F, sep="\t", fill=T, check.names=F, quote="") #318
length(as.character(og_data$V1)) #318

oncogene_list <- unique(c(as.character(pathway_genes[which(pathway_genes$`OG/TSG` == "OG"), 1]), as.character(og_data$V1))) #392

og_gene_DEG <- c()
for(i in 1:nrow(DEG)){
  if(DEG[i,1] %in% oncogene_list){
    og_gene_DEG <- rbind(og_gene_DEG, DEG[i,])
  }
}

write.table( og_gene_DEG, "og_gene_DEG_padjust.txt", sep="\t", row.names=T, col.names=T, quote=F )

#volcano plot
sub_DEG <- DEG[na.omit(match(oncogene_list, DEG$symbol)), ]
plot(log2(sub_DEG[,"fc_IDHn"]), -log10(sub_DEG[,"padjust_IDHn"]), 
     xlab="log2(fold change)", ylab="-log10(p.adjust)", pch=20, col="gray", cex.lab=1.3, 
     xlim=c(-6, 6), 
     las=1, bty="l",main = "All")
abline(h=-log10(0.05), lty=2)
abline(v=log2(0.5), lty=2)
abline(v=log2(2), lty=2)
#79
sub_up <- sub_DEG[ which(sub_DEG[,"padjust_IDHn"] <= 0.05 & sub_DEG[,"fc_IDHn"] >= 2),]
#21
sub_down <- sub_DEG[ which(sub_DEG[,"padjust_IDHn"] <= 0.05 & sub_DEG[,"fc_IDHn"] <= 0.5),]

phUMR_OGs <- c("MAML2", "ZEB1", "MYCL", "BCL6", "PDGFRA", "CCND1", "SYK", "MAP3K1", "CXCR4", "RUNX1", "LIMD1", "SIX1", "AJUBA", "HOXA13", "ACKR3", "TET2", "MSI2", "TET1", "ETV6")
phUMR_OG_DEG <- DEG[na.omit(match(phUMR_OGs, DEG$symbol)), ]
points(log2(phUMR_OG_DEG[,"fc_IDHn"]), -log10(phUMR_OG_DEG[,"padjust_IDHn"]), pch=20, col=red)

pdf_file="E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\pdf\\oncogene_phUMRs.pdf"
dev.print(pdf, file=pdf_file, width=4, height=4)

partial_down <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\genes\\partial_down_genes.bed",
                  header=F, sep="\t", fill=T, check.names=F, quote="")
full_down <- read.table("E:\\IDH_Hyper\\Results\\4_Gene_Expression\\202210_genes_GO_volcano\\genes\\full_down_genes.bed",
                           header=F, sep="\t", fill=T, check.names=F, quote="")

phUMR_down_DEG <- sub_DEG[na.omit(match(partial_down$V5, sub_DEG$symbol)), ]
fhUMR_down_DEG <- sub_DEG[na.omit(match(full_down$V5, sub_DEG$symbol)), ]


points(log2(phUMR_OG_DEG[,"fc_IDHn"]), -log10(phUMR_OG_DEG[,"padjust_IDHn"]), pch=20, col=red)

#TSG
tsg_data <- read.table("E:\\IDH_Hyper\\Data\\COSMIC\\COSMIC_TSG.txt",
                       header=F, sep="\t", fill=T, check.names=F, quote="") #316
#386
tsg_list <- unique(c(as.character(pathway_genes[which(pathway_genes$`OG/TSG` == "TSG"), 1]), as.character(tsg_data$V1))) 

tsg_gene_DEG <- c()
for(i in 1:nrow(DEG)){
  if(DEG[i,1] %in% tsg_list){
    tsg_gene_DEG <- rbind(tsg_gene_DEG, DEG[i,])
  }
}

write.table( tsg_gene_DEG, "tsg_gene_DEG_padjust.txt", sep="\t", row.names=T, col.names=T, quote=F )

rownames(count) <- c("partial_up", "partial_down", "full_up", "full_down")
colnames(count) <- c("OG", "TSG")

og_pvalue <- c()
for(i in 1:nrow(count)){
  og_pvalue <- c(og_pvalue, phyper(count[i,1], 244, 56202-244, 278, lower.tail=F)) #242
}
tsg_pvalue <- c()
for(i in 1:nrow(count)){
  tsg_pvalue <- c(tsg_pvalue, phyper(count[i,2], 242, 56202-242, 124, lower.tail=F)) #242
}

pvalue <- cbind(og_pvalue, tsg_pvalue)
rownames(pvalue) <- c("partial_up", "partial_down", "full_up", "full_down")
colnames(pvalue) <- c("OG", "TSG")

pvalue_log <- -log10(pvalue)

library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggprism)

gg_count <- melt(pvalue_log)
colnames(gg_count) <- c("hyper_type", "Gene_type", "value")

ggplot(gg_count, aes(Gene_type, value, color=Gene_type, fill=Gene_type)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_grid(~hyper_type,scales = 'free') +
  labs(x = "", y = "-log(p-value)") +
  theme_prism() 

dev.print(pdf, file="OG_TSG_pvalue.pdf", width=6, height=4)
