
blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH GBM
red = rgb(255, 45, 41, max=255) #WT GBM

###1. Normal - 150 -- UMR -- 组蛋白修饰
samples <- c("149_H3K27ac", "149_H3K27me3", "149_H3K36me3", "149_H3K4me1", "149_H3K4me3", "149_H3K9me3",
             "150_H3K27ac", "150_H3K27me3", "150_H3K36me3", "150_H3K4me1", "150_H3K4me3", "150_H3K9me3")

for(i in samples){
  hist_distribution(i)
}

#partial_K27me3_bed <- read.table("signal/partial_150_H3K27me3.bed", header=F, sep="\t", fill=T, check.names=F, quote="")
#plot(partial_K27me3_bed$V3 - partial_K27me3_bed$V2, partial_K27me3_bed$V5)
#full_K27me3_bed <- read.table("signal/full_150_H3K27me3.bed", header=F, sep="\t", fill=T, check.names=F, quote="")
#plot(full_K27me3_bed$V3 - full_K27me3_bed$V2, full_K27me3_bed$V5)

hist_distribution_fig("150_H3K4me3", "f_p", -5, 100)
hist_distribution_fig("150_H3K27ac", "f_p", -5, 100)
hist_distribution_fig("150_H3K4me1", "p_f", -2, 30)
hist_distribution_fig("150_H3K9me3", "p_f", -2, 30)
hist_distribution_fig("150_H3K27me3", "f_p", -2, 30)
hist_distribution_fig("150_H3K36me3", "p_f", -2, 30)

hist_distribution_fig <- function(signal_id , type, xmin, xmax){ #type: p_f, f_p
  full_file <- paste("signal/full_", signal_id, ".bed", sep="")
  partial_file <- paste("signal/partial_", signal_id, ".bed", sep="")
  
  full_bed <- read.table(full_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  partial_bed <- read.table(partial_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  
  if(type == "p_f"){
    plot(density( partial_bed$V5 ) ,
         xlab = "signal", ylab = "Density",
         xlim =c(xmin, xmax),
         col=red, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
    lines(density(full_bed$V5), ylab="", xlab="", col=blue, lwd=2, cex.lab=1.3, las=1)
    dev.print(pdf, file=paste("Rplot/", signal_id, ".pdf", sep=""), width=4, height=6)
  }
  if(type == "f_p"){
    plot(density( full_bed$V5 ) ,
         xlab = "signal", ylab = "Density",
         xlim =c(xmin, xmax),
         col=blue, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
    lines(density(partial_bed$V5), ylab="", xlab="", col=red, lwd=2, cex.lab=1.3, las=1)
    dev.print(pdf, file=paste("Rplot/", signal_id, ".pdf", sep=""), width=4, height=6)
  }
  
  return(0)
}


hist_distribution <- function(signal_id){ #注：matrix 不能有NA
  full_file <- paste("signal/full_", signal_id, ".bed", sep="")
  partial_file <- paste("signal/partial_", signal_id, ".bed", sep="")
  
  full_bed <- read.table(full_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  partial_bed <- read.table(partial_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  
  plot(density( partial_bed$V5 ) ,
       xlab = "signal", ylab = "Density",
       col=red, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
  lines(density(full_bed$V5), ylab="", xlab="", col=blue, lwd=2, cex.lab=1.3, las=1)
  dev.print(pdf, file=paste("Rplot/", signal_id, ".pdf", sep=""), width=4, height=6)
  return(0)
}

###2. IDH - 组蛋白修饰
hist_dist_IDH_fig("AK076_H3K4me3", "p_f", -5, 100)
hist_dist_IDH_fig("AK076_H3K27ac", "p_f", -5, 100)
hist_dist_IDH_fig("AK076_H3K4me1", "p_f", -2, 30)
hist_dist_IDH_fig("AK076_H3K9me3", "p_f", -2, 30)
hist_dist_IDH_fig("AK076_H3K27me3", "p_f", -2, 30)
hist_dist_IDH_fig("AK076_H3K36me3", "p_f", -2, 30)

hist_dist_IDH_fig("AK124_H3K4me3", "p_f", -5, 100)
hist_dist_IDH_fig("AK124_H3K27ac", "p_f", -5, 100)
hist_dist_IDH_fig("AK124_H3K4me1", "p_f", -2, 30)
hist_dist_IDH_fig("AK124_H3K9me3", "p_f", -2, 30)
hist_dist_IDH_fig("AK124_H3K27me3", "p_f", -2, 30)
hist_dist_IDH_fig("AK124_H3K36me3", "p_f", -2, 30)

hist_dist_IDH_fig <- function(signal_id , type, xmin, xmax){ #type: p_f, f_p
  full_file <- paste("signal_IDH/full_hyper_", signal_id, ".bed", sep="")
  partial_file <- paste("signal_IDH/partial_hyper_", signal_id, ".bed", sep="")
  
  full_bed <- read.table(full_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  partial_bed <- read.table(partial_file, header=F, sep="\t", fill=T, check.names=F, quote="")
  
  if(type == "p_f"){
    plot(density( partial_bed$V5 ) ,
         xlab = "signal", ylab = "Density",
         xlim =c(xmin, xmax),
         col=red, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
    lines(density(full_bed$V5), ylab="", xlab="", col=blue, lwd=2, cex.lab=1.3, las=1)
    dev.print(pdf, file=paste("Rplot_IDH/", signal_id, ".pdf", sep=""), width=4, height=6)
  }
  if(type == "f_p"){
    plot(density( full_bed$V5 ) ,
         xlab = "signal", ylab = "Density",
         xlim =c(xmin, xmax),
         col=blue, lwd=2, cex.lab=1.3, cex.axis=1, las=1, main="", bty="l", font.axis=2)
    lines(density(partial_bed$V5), ylab="", xlab="", col=red, lwd=2, cex.lab=1.3, las=1)
    dev.print(pdf, file=paste("Rplot_IDH/", signal_id, ".pdf", sep=""), width=4, height=6)
  }
  
  return(0)
}

