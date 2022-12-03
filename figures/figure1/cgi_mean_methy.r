#boxplot of CGI mean methylation for normal and IDH glioma
blue = rgb(0, 140, 218, max=255) #brain normal
orange = rgb(243, 146, 0, max=255) #IDH glioma

brain_CGI_methy <- read.table("methy_matrix\\CGI_brain_normal_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")
gbm_CGI_methy <- read.table("methy_matrix\\CGI_gbm_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")

brain_IDH_Hyper_methy <- read.table("methy_matrix\\IDH_Hyper_brain_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")
gbm_IDH_Hyper_methy <- read.table("methy_matrix\\IDH_Hyper_GBM_methy_matrix.txt", header=T, sep="\t", fill=T, check.names=F, quote="")

#gbm_CGI_methy[1:5, 1:5]
#colnames(gbm_CGI_methy)[c(8, 11, 17, 18, 22, 25, 32, 33, 36, 54, 56, 62, 64, 65, 66)] #IDH samples
colnames(gbm_CGI_methy)[c(8, 11, 17, 18, 22, 25, 32, 33, 36, 54, 56, 62, 64, 65, 66)] <- "IDH"
#check - right
normal_index <- c(4:ncol(brain_CGI_methy))
IDH_index <- c(8, 11, 17, 18, 22, 25, 32, 33, 36, 54, 56, 62, 64, 65, 66)

#boxplot
#index - line - 1
CGI_index <- 19914 #AFF1 CGI, chr4,87855873,87857204
IDH_Hyper_index <- 3276 #AFF1 IDH Hyper, chr4	87856882	87857184
CGI_index <- 25508 #MLLT3 CGI, chr9,20620728,20621862
IDH_Hyper_index <- 4244 #MLLT3 IDH Hyper, chr9	20619900	20620425

CGI_index <- 1459 #SYT6 CGI, chr1	114695136	114696672
IDH_Hyper_index <- 254 #SYT6 IDH Hyper, chr1	114696191	114697369

CGI_index <- 9402 #TPPP3 CGI, chr16	67427284	67428950
IDH_Hyper_index <- 1619 #TPPP3 IDH Hyper, chr16	67427048	67427843

CGI_index <- 4502 #CCND1 CGI, chr11	69451136	69458596
IDH_Hyper_index <- 798 #CCND1 IDH Hyper, chr11	69451557	69452302

brain_CGI_methy[CGI_index,c(1:3)]
brain_IDH_Hyper_methy[IDH_Hyper_index,c(1:3)]

brain_methy <- na.omit(as.numeric(brain_CGI_methy[CGI_index, normal_index]))
IDH_methy <- na.omit(as.numeric(gbm_CGI_methy[CGI_index, IDH_index]))

adjborder_brain_methy <- na.omit(as.numeric(brain_IDH_Hyper_methy[IDH_Hyper_index, normal_index]))
adjborder_IDH_methy <- na.omit(as.numeric(gbm_IDH_Hyper_methy[IDH_Hyper_index, IDH_index]))

boxplot(brain_methy, IDH_methy, adjborder_brain_methy, adjborder_IDH_methy,
        col=c(blue, red), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 2, 
        ylim=c(0, 1), names=c("Normal Brain", "IDH Glioma", "Normal Brain", "IDH Glioma"),
        main="CCND1 Methylation changes")
mean(IDH_methy) - mean(brain_methy) #0.1429333
mean(adjborder_IDH_methy) - mean(adjborder_brain_methy) #0.6410857
dev.print(pdf, file="CCND1_methy.pdf", width=8, height=8)

boxplot(brain_methy, IDH_methy, WT_methy, adjborder_brain_methy, adjborder_IDH_methy, adjborder_WT_methy,
        col=c(blue, red, orange), las=1, bty="l",
        cex.axis = 1.5, cex.lab=1.5, cex.main = 2, 
        ylim=c(0, 1), names=c("Normal Brain", "IDH Glioma", "WT Glioma", "Normal Brain", "IDH Glioma", "WT Glioma"),
        main="TPPP3 Methylation changes")


###comparison
dim(na.omit(brain_CGI_methy))
CGI_NA_counts <- c()
for(i in 1:nrow(brain_CGI_methy)){
  na_stat <- is.na(as.numeric(brain_CGI_methy[i,normal_index]))
  CGI_NA_counts <- c( CGI_NA_counts, length( which(na_stat == FALSE) ) )
}
length( which(CGI_NA_counts > 45) ) #more than 45 samples methylation

overlaped_CGI <- read.table("methy_matrix\\overlaped_CGI.bed", header=F, sep="\t", fill=T, check.names=F, quote="")

all_cgi_names <- c()
for(i in 1:nrow(brain_CGI_methy)){
  all_cgi_names <- c(all_cgi_names, paste(brain_CGI_methy[i,1:3],collapse=","))
}
rownames(brain_CGI_methy) <- all_cgi_names

overlaped_cgi_names <- c()
for(i in 1:nrow(overlaped_CGI)){
  overlaped_cgi_names <- c(overlaped_cgi_names, paste(overlaped_CGI[i,1:3],collapse=","))
}

dim(brain_CGI_methy[overlaped_cgi_names,])

length(setdiff(all_cgi_names, overlaped_cgi_names))

only_CGI_brain_methy <- brain_CGI_methy[setdiff(all_cgi_names, overlaped_cgi_names), ]
only_CGI_brain_methy_arr <- c()
for(i in 1:nrow(only_CGI_brain_methy)){
  only_CGI_brain_methy_arr <- c(only_CGI_brain_methy_arr, as.numeric(only_CGI_brain_methy[i,normal_index]))
}
plot(density(na.omit(only_CGI_brain_methy_arr))) #binormal

only_CGI_brain_methy_range <- c()
for(i in 1:nrow(only_CGI_brain_methy)){
  methy <- na.omit(as.numeric(only_CGI_brain_methy[i,normal_index]))
  only_CGI_brain_methy_range <- c(only_CGI_brain_methy_range, max(methy) - min(methy))
}
plot(density(na.omit(only_CGI_brain_methy_range))) #binormal

dim(only_CGI_brain_methy [ which(only_CGI_brain_methy_range <= 0.3), normal_index])

only_CGI_brain_methy_mean

library(pheatmap)
pheatmap( na.omit(only_CGI_brain_methy [ which(only_CGI_brain_methy_mean >= 0.6), normal_index]),
          color = colorRampPalette(c("blue", "white", "red"))(100),
          show_rownames = T, show_colnames = T,
          cluster_rows = T, cluster_cols = T,
          clustering_method = "ward.D2")

only_CGI_brain_methy2 <- only_CGI_brain_methy [ which(only_CGI_brain_methy_mean >= 0.6), normal_index]
only_CGI_brain_methy_arr2 <- c()
for(i in 1:nrow(only_CGI_brain_methy2)){
  only_CGI_brain_methy_arr2 <- c(only_CGI_brain_methy_arr2, as.numeric(only_CGI_brain_methy[i,normal_index]))
}
plot(density(na.omit(only_CGI_brain_methy_arr2))) 

