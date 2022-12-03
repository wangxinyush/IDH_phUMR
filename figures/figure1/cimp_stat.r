#compared with CIMP
blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH GBM
red = rgb(255, 45, 41, max=255) #WT GBM

brain_CGI_methy <- read.table("Data\\CGI_brain_normal_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")
gbm_CGI_methy <- read.table("Data\\CGI_gbm_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")

colnames(gbm_CGI_methy) <- c("chr", "start", "end", rep("WT", 66))
colnames(gbm_CGI_methy)[c(8, 11, 17, 18, 22, 25, 32, 33, 36, 54, 56, 62, 64, 65, 66)] <- "IDH"
#check - right
normal_index <- c(4:ncol(brain_CGI_methy))
IDH_index <- c(8, 11, 17, 18, 22, 25, 32, 33, 36, 54, 56, 62, 64, 65, 66)
WT_index <- which(colnames(gbm_CGI_methy) == "WT")

rownames(brain_CGI_methy) <- paste(gbm_CGI_methy[,1],gbm_CGI_methy[,2],gbm_CGI_methy[,3], sep=",")
rownames(gbm_CGI_methy) <- paste(gbm_CGI_methy[,1],gbm_CGI_methy[,2],gbm_CGI_methy[,3], sep=",")

#partial_CGI / full_CGI
partial_CGI <- read.table("Data\\partial_CGI.bed", header=F, sep="\t", fill=T, check.names=F, quote="")
full_CGI <- read.table("Data\\full_CGI.bed", header=F, sep="\t", fill=T, check.names=F, quote="")
partial_CGI_str <- paste(partial_CGI[,1],partial_CGI[,2],partial_CGI[,3], sep=",")
full_CGI_str <- paste(full_CGI[,1],full_CGI[,2],full_CGI[,3], sep=",")

#CIMP
partial_CIMP_stat <- c()
for (i in 1:length(partial_CGI_str)){
  CGI_str <- partial_CGI_str[i]
  
  IDH_methy <- na.omit(as.numeric(gbm_CGI_methy[CGI_str, IDH_index]))
  WT_methy <- na.omit(as.numeric(gbm_CGI_methy[CGI_str, WT_index]))
  norm_methy <- na.omit(as.numeric(brain_CGI_methy[CGI_str, normal_index]))
  
  if( length(IDH_methy) >= 3 && length(WT_methy) >= 3 ){
    IDH_WT_p <- t.test(IDH_methy, WT_methy)$p.value
  }
  else{
    IDH_WT_p <- 1
  }
  if( length(IDH_methy) >= 3 && length(norm_methy) >= 3 ){
    IDH_Norm_p <- t.test(IDH_methy, norm_methy)$p.value
  }
  else{
    IDH_Norm_p <- 1
  }
  
  partial_CIMP_stat <- rbind(partial_CIMP_stat, c(mean(IDH_methy), mean(WT_methy), mean(norm_methy), IDH_WT_p, IDH_Norm_p))
}
colnames(partial_CIMP_stat) <- c("IDH_mean", "WT_mean", "Normal_mean", "IDH_WT_p", "IDH_Norm_p")
IDH_WT_p_adjust <- p.adjust(partial_CIMP_stat[,4], method = "BH", n = nrow(partial_CIMP_stat))
IDH_Norm_p_adjust <- p.adjust(partial_CIMP_stat[,5], method = "BH", n = nrow(partial_CIMP_stat))
partial_CIMP_stat <- cbind(partial_CIMP_stat, IDH_WT_p_adjust, IDH_Norm_p_adjust)

#full CIMP
full_CIMP_stat <- c()
for (i in 1:length(full_CGI_str)){
  CGI_str <- full_CGI_str[i]
  
  IDH_methy <- na.omit(as.numeric(gbm_CGI_methy[CGI_str,IDH_index]))
  WT_methy <- na.omit(as.numeric(gbm_CGI_methy[CGI_str,WT_index]))
  norm_methy <- na.omit(as.numeric(brain_CGI_methy[CGI_str,normal_index]))
  
  if( length(IDH_methy) >= 3 && length(WT_methy) >= 3 ){
    IDH_WT_p <- t.test(IDH_methy, WT_methy)$p.value
  }
  else{
    IDH_WT_p <- 1
  }
  if( length(IDH_methy) >= 3 && length(norm_methy) >= 3 ){
    IDH_Norm_p <- t.test(IDH_methy, norm_methy)$p.value
  }
  else{
    IDH_Norm_p <- 1
  }
  
  full_CIMP_stat <- rbind(full_CIMP_stat, c(mean(IDH_methy), mean(WT_methy), mean(norm_methy), IDH_WT_p, IDH_Norm_p))
}
colnames(full_CIMP_stat) <- c("IDH_mean", "WT_mean", "Normal_mean", "IDH_WT_p", "IDH_Norm_p")
IDH_WT_p_adjust <- p.adjust(full_CIMP_stat[,4], method = "BH", n = nrow(full_CIMP_stat))
IDH_Norm_p_adjust <- p.adjust(full_CIMP_stat[,5], method = "BH", n = nrow(full_CIMP_stat))
full_CIMP_stat <- cbind(full_CIMP_stat, IDH_WT_p_adjust, IDH_Norm_p_adjust)

#stat
#CIMP
#380
partial_CIMP_len <- length(which(partial_CIMP_stat[,"IDH_mean"] - partial_CIMP_stat[,"WT_mean"] >= 0.2 
             & partial_CIMP_stat[, "IDH_WT_p_adjust"] < 0.05) )
#658
full_CIMP_len <- length(which(full_CIMP_stat[,"IDH_mean"] - full_CIMP_stat[,"WT_mean"] >= 0.2
             & full_CIMP_stat[, "IDH_WT_p_adjust"] < 0.05) )

#346
length(which(partial_CIMP_stat[,"IDH_mean"] - partial_CIMP_stat[,"WT_mean"] >= 0.2 
             & partial_CIMP_stat[, "IDH_WT_p_adjust"] < 0.05
             & partial_CIMP_stat[,"IDH_mean"] - partial_CIMP_stat[,"Normal_mean"] >= 0.2
             & partial_CIMP_stat[, "IDH_Norm_p_adjust"] < 0.05 ) )
#653
length(which(full_CIMP_stat[,"IDH_mean"] - full_CIMP_stat[,"WT_mean"] >= 0.2
             & full_CIMP_stat[, "IDH_WT_p_adjust"] < 0.05
             & full_CIMP_stat[,"IDH_mean"] - full_CIMP_stat[,"Normal_mean"] >= 0.2
             & full_CIMP_stat[, "IDH_Norm_p_adjust"] < 0.05 ) )

percentage <- 100*c(1 - partial_CIMP_len/length(partial_CGI_str), 1 - full_CIMP_len/length(full_CGI_str))
names(percentage) <- c("Partially", "Fully")
barplot(percentage, 
        #ylim = c(0, 100), 
        beside=T, col=c(red, "gray"), border=T, las=1, cex.lab=1.3, ylab="Percentage")

