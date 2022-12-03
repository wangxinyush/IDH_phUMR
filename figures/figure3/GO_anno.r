#GO annotation of phUMR, fhUMR related genes

blue = rgb(15, 112, 183, max=255)
red = rgb(190, 21, 34, max=255)
orange = rgb(243, 146, 0, max=255)

##Data input - partially/fully hyper genes
promoter_partial_up <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\DEG_regions\\promoter_partial_up_anno.txt",
                             header=F, sep="\t", fill=T, check.names=F, quote="")
promoter_partial_down <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\DEG_regions\\promoter_partial_down_anno.txt",
                         header=F, sep="\t", fill=T, check.names=F, quote="")
promoter_full_down <- read.table("E:\\IDH_Hyper\\Results\\2_Annotation_IDH_Hyper\\regions\\DEG_regions\\promoter_full_down_anno.txt",
                                    header=F, sep="\t", fill=T, check.names=F, quote="")

partial_up_genes <- unique( as.character(promoter_partial_up$V6) )
partial_down_genes <- unique( as.character(promoter_partial_down$V6) )
full_down_genes <- unique( as.character(promoter_full_down$V6) )

write.table( partial_up_genes, "genes\\partial_up_genes_symbol.txt", sep="\t", row.names=F, col.names=F, quote=F )
write.table( partial_down_genes, "genes\\partial_down_genes_symbol.txt", sep="\t", row.names=F, col.names=F, quote=F )
write.table( full_down_genes, "genes\\full_down_genes_symbol.txt", sep="\t", row.names=F, col.names=F, quote=F )

#length(intersect(partial_down_genes, full_down_genes)) #0

partial_up_genes_entrez <- symbol_to_entrez(partial_up_genes)
partial_down_genes_entrez <- symbol_to_entrez(partial_down_genes)
full_down_genes_entrez <- symbol_to_entrez(full_down_genes)

#GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(annotate)

#1) GO
partial_up_ego <- enrichGO(gene=partial_up_genes_entrez, 'org.Hs.eg.db', keyType= "ENTREZID",ont="BP", pvalueCutoff=1, qvalueCutoff=1, minGSSize = 0, maxGSSize = 5000)
partial_down_ego <- enrichGO(gene=partial_down_genes_entrez, 'org.Hs.eg.db', keyType= "ENTREZID",ont="BP", pvalueCutoff=1, qvalueCutoff=1, minGSSize = 0, maxGSSize = 5000)
full_down_ego <- enrichGO(gene=full_down_genes_entrez, 'org.Hs.eg.db', keyType= "ENTREZID",ont="BP", pvalueCutoff=1, qvalueCutoff=1, minGSSize = 0, maxGSSize = 5000)

write.table( partial_up_ego@result, "ego_results\\promoter_partial_up_ego_results.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( partial_down_ego@result, "ego_results\\promoter_partial_down_ego_results.txt", sep="\t", row.names=F, col.names=T, quote=F )
write.table( full_down_ego@result, "ego_results\\promoter_full_down_ego_results.txt", sep="\t", row.names=F, col.names=T, quote=F )

partial_up_all_table <- partial_up_ego@result
partial_down_all_table <- partial_down_ego@result
full_down_all_table <- full_down_ego@result

plot_GO_ID <- c("GO:0045595", "GO:0048589", "GO:0007423", "GO:0045664", "GO:0006811",
                "GO:0007268", "GO:0042391", "GO:0007399", "GO:0034762", "GO:0006775", "GO:1902667", "GO:0099612") 

dim(partial_up_all_table)
dim(partial_down_all_table)
dim(full_down_all_table)
test_p <- head(  full_down_all_table$pvalue  )
-log10(p.adjust(test_p, method = "BH", n=200))

plot_GO_ID <- rev(plot_GO_ID)
partial_full_up_down_GO <- rbind(partial_up_all_table[plot_GO_ID, "p.adjust"], 
                                 partial_down_all_table[plot_GO_ID, "p.adjust"],
                                 full_down_all_table[plot_GO_ID, "p.adjust"])
colnames(partial_full_up_down_GO) <- partial_up_all_table[plot_GO_ID, "Description"]
rownames(partial_full_up_down_GO) <- c("Partially_up", "Partially_down", "Fully_down")

barplot( -log10(partial_full_up_down_GO), beside = T, horiz = T,
         xlim = c(0, 10),
         col=c(red, blue, "gray"), names.arg=colnames(partial_full_up_down_GO), las=2, xlab = "-log10(adjusted p-value)")

abline(v=-log10(0.05), lty=2)
dev.print(pdf, file="partial_full_up_down_GO.pdf", width=10, height=8)

#2) KEGG
#use Metascape
partial_up_KEGG <- read.table("Metascape\\partial_up_KEGG.txt",
                             header=T, sep="\t", fill=T, check.names=F, quote="")
partial_down_KEGG <- read.table("Metascape\\partial_down_KEGG.txt",
                              header=T, sep="\t", fill=T, check.names=F, quote="")
full_down_KEGG <- read.table("Metascape\\full_down_KEGG.txt",
                              header=T, sep="\t", fill=T, check.names=F, quote="")

#partial_up_KEGG_padjust <- p.adjust(10^partial_up_KEGG$LogP, method = "BH", n=52)
#partial_down_KEGG_padjust <- p.adjust(10^partial_down_KEGG$LogP, method = "BH", n=52)
#full_down_KEGG_padjust <- p.adjust(10^full_down_KEGG$LogP, method = "BH", n=52)
#partial_up_KEGG <- cbind(partial_up_KEGG, LogPadjust=-log10(partial_up_KEGG_padjust))
#intersect(partial_up_KEGG$`Term_ID`, partial_down_KEGG$`Term_ID`) #"hsa04010" "hsa04360"
#intersect(partial_down_KEGG$`Term_ID`, full_down_KEGG$`Term_ID`) #null

pathway_name <- c("Hippo signaling pathway", "Pathways in cancer", "Wnt signaling pathway", 
                  "TGF-beta signaling pathway", "Axon guidance", "Apoptosis", "Morphine addiction", 
                  "GABAergic synapse", "Neuroactive ligand-receptor interaction", "Ras signaling pathway", 
                  "Dopaminergic synapse", "PPAR signaling pathway")
pathway_name <- rev(pathway_name)
pathway_padjust <- c()
for(i in 1:length(pathway_name)){
  pathway <- pathway_name[i]
  
  partial_up_count_str <- partial_up_KEGG[ which(partial_up_KEGG$ID == pathway)[1] , "InTerm_InList"]
  partial_up_pathway_count <- as.numeric( strsplit(as.character(partial_up_count_str), "/")[[1]][1] )
  pathway_gene_count1 <- as.numeric( strsplit(as.character(partial_up_count_str), "/")[[1]][2] )
  
  partial_down_count_str <- partial_down_KEGG[ which(partial_down_KEGG$ID == pathway)[1] , "InTerm_InList"]
  partial_down_pathway_count <- as.numeric( strsplit(as.character(partial_down_count_str), "/")[[1]][1] )
  pathway_gene_count2 <- as.numeric( strsplit(as.character(partial_down_count_str), "/")[[1]][2] )
  
  full_down_count_str <- full_down_KEGG[ which(full_down_KEGG$ID == pathway)[1] , "InTerm_InList"]
  full_down_pathway_count <- as.numeric( strsplit(as.character(full_down_count_str), "/")[[1]][1] )
  pathway_gene_count3 <- as.numeric( strsplit(as.character(full_down_count_str), "/")[[1]][2] )
  
  pathway_gene_count <- max(na.omit(c(pathway_gene_count1, pathway_gene_count2, pathway_gene_count3)))
  if( is.na(partial_up_pathway_count) ){
    partial_up_pathway_count <- 0
  }
  partial_up_p <- phyper(partial_up_pathway_count, pathway_gene_count,35044-pathway_gene_count,130,lower.tail=F) #
  
  if( is.na(partial_down_pathway_count) ){
    partial_down_pathway_count <- 0
  }
  partial_down_p <- phyper(partial_down_pathway_count, pathway_gene_count,35044-pathway_gene_count,236,lower.tail=F) #
  
  if( is.na(full_down_pathway_count) ){
    full_down_pathway_count <- 0
  }
  full_down_p <- phyper(full_down_pathway_count, pathway_gene_count,35044-pathway_gene_count,66,lower.tail=F) #
  
  padjust <- p.adjust(c(partial_up_p, partial_down_p, full_down_p), method = "BH")
  pathway_padjust <- rbind(pathway_padjust, padjust)
}
rownames(pathway_padjust) <- pathway_name
colnames(pathway_padjust) <- c("partial_up", "partial_down", "full_down")

#pathway_padjust["PPAR signaling pathway", "full_down"] <- full_down_KEGG

barplot( t(-log10(pathway_padjust)), beside = T, horiz = T, 
         xlim = c(0,14),
         col=c(red, blue, "gray"),
         names.arg=rownames(pathway_padjust),
        las=2, xlab = "-log10(adjusted p-value)")
abline(v=-log10(0.05), lty=2)
dev.print(pdf, file="partial_full_up_down_KEGG.pdf", width=10, height=8)

#clusterprofiler ID conservation
get_symbol_gene <- function(peakAnnoList_anno, dis_cutoff){ #peakAnnoList$`WT_WT`, 10000
  require(org.Hs.eg.db)
  require(annotate)
  anno <- as.data.frame(peakAnnoList_anno)
  gene_id <- anno$geneId[which( abs(anno$distanceToTSS ) <= dis_cutoff)]
  return(getSYMBOL(as.character(gene_id), data='org.Mm.eg'))
}
get_entrez_gene <- function(peakAnnoList_anno, dis_cutoff){ #peakAnnoList$`WT_WT`, 10000
  require(org.Hs.eg.db)
  require(annotate)
  anno <- as.data.frame(peakAnnoList_anno)
  gene_id <- anno$geneId[which( abs(anno$distanceToTSS ) <= dis_cutoff)]
  return(gene_id)
}
symbol_to_entrez <- function(gene_symbol){
  library(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  result <- select(hs, keys = gene_symbol, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
  return(result$ENTREZID)
}

