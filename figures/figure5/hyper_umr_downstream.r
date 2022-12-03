#This program was used to plot delta histone modifications.

blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH GBM
red = rgb(255, 45, 41, max=255) #WT GBM
purple = rgb(153, 0, 153, max=255) #
green = rgb(12, 99, 52, max=255) #CGI

cgi_blue = rgb(57, 127, 185, max=255) #CGI
shore_purple = rgb(152, 79, 159, max=255) #shore
shelf_green = rgb(80, 175, 73, max=255) #shelf
sea_gray = rgb(137, 137, 137, max=255) #sea

onco_histone_plot("b149", "g076", -3.5, 3.5) #0.5
onco_histone_plot("b150", "g213", -3.5, 3.5)

onco_histone_plot<- function(sample1, sample2, ymin, ymax){ #注：matrix 不能有NA
  # sample1 <- "b149"
  # sample2 <- "g076"
  # ymin <- -5
  # ymax <- 5
  
  histone <- c("H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K27me3", "H3K9me3")
  region <- c("hyper", "umr", "downstream")
  
  fc_matrix <- c()
  for(j in 1:length(histone)){
  for(i in 1:length(region)){
    b_signal_file <- paste("histone/", sample1, "_", histone[j], "_", region[i], ".bed", sep="")
    g_signal_file <- paste("histone/", sample2, "_", histone[j], "_", region[i], ".bed", sep="")
    b_signal_data <- read.table(b_signal_file, header=F, sep="\t", fill=T, check.names=F, quote="")
    g_signal_data <- read.table(g_signal_file, header=F, sep="\t", fill=T, check.names=F, quote="")
    
    b_signal <- b_signal_data$V5
    g_signal <- g_signal_data$V5
    
    b_signal[ which(b_signal <= 1) ] <- 0.5
    g_signal[ which(g_signal <= 1) ] <- 0.5
    
    fc <- log2(g_signal / b_signal)
    fc_matrix <- cbind(fc_matrix, fc)
  }
  }
  
  boxplot(fc_matrix[,1], fc_matrix[,2], fc_matrix[,3],#hyper, umr, downstream
          fc_matrix[,4], fc_matrix[,5], fc_matrix[,6],
          fc_matrix[,7], fc_matrix[,8], fc_matrix[,9],
          fc_matrix[,10], fc_matrix[,11], fc_matrix[,12],
          fc_matrix[,13], fc_matrix[,14], fc_matrix[,15],
          fc_matrix[,16], fc_matrix[,17], fc_matrix[,18],
          col=c(red, cgi_blue, "gray"), las=1, bty="l", outline = F, ylim=c(ymin, ymax),
          ylab="log2(FC)", 
          #names = c("Hyper", "UMR", "Downstream"),
          cex.axis = 1.5, cex.lab=1.5, cex.main = 2 )
  dev.print(pdf, file=paste("pdf/",sample1,"vs",sample2, ymin,"_histone.pdf",sep=""), width=10, height=6)
  return(0)
}
