#This program was used to identify three methylation states sites by using HMM.

##Use bed_mean_methy on the Server to map the methylation positions of 9 files(brain3gbm6) to the reference CGI / reference UMR
##bed_mean_methy -r /share/pub/wangxy/Annotation/CGI/hg19_CGI.bed -w brain3gbm6.txt -o cgibedmean.bed
cgi_bedmean <- read.table("E:\\IDH\\R\\cgibedmean.bed")
####Determine the states of each CGI
#dim(cgi_bedmean)  #27719  12
cgi_mean <- na.omit(cgi_bedmean) #remove NA lines
#dim(cgi_mean)  #25614  12
cgi_mean_ma <- as.matrix(cgi_mean)[-1,]  #convert data to matrix
cgi_states <- matrix(nrow = nrow(cgi_mean_ma),ncol = (ncol(cgi_mean_ma)+1))
hyper <- c() #Prepare to generate sets of three states
nodiff <- c()
hypo <- c()
for (i in 1:nrow(cgi_mean_ma)){
  a <- as.numeric(cgi_mean_ma[i,10:12]) #set a as idhmut methylation
  b <- as.numeric(cgi_mean_ma[i,4:6]) #set b as normal methylation
  mean <-mean(a)-mean(b)  #Calculate the respective mean
  c <- which(a==a[1])   
  d <- which(b==b[1])
  if((length(c)!=3) & (length(d)!=3)){ #remove constant data
    p <- t.test(a,b)$p.value  
    if((mean >= 0.2)&(p<0.05)){
      hyper <- as.numeric(c(hyper,mean))
    }
    else if((mean <= -0.2)&(p<0.05)){
      hypo <- as.numeric(c(hypo,mean))
    }else {
      nodiff <- as.numeric(c(nodiff,mean))
    }
  }
}

#length(hyper) #3082
#length(nodiff)  # 20496
#length(hypo)  # 8

####Emission probability matrix
options(digits = 3)
means <- c(mean(hyper),mean(nodiff),mean(hypo))
means #0.4258  0.0443 -0.2675
sds <- c(sd(hyper),sd(nodiff),sd(hypo))
sds  #0.1401 0.0973 0.0414

####transition probability matrix
###Hyper→Hyper,no-diff→no-diff,Hypo→Hypo
##bed_mean_methy -r ref_UM.bed -w brain3gbm6.txt -o umrbedmean.bed
umr_bedmean <- read.table("E:\\IDH\\R\\umrbedmean.bed")
#dim(umr_bedmean)  #26065    12
umr_mean <- na.omit(umr_bedmean) #remove NA lines
#dim(umr_mean)  #25008  12
umr_mean_ma <- as.matrix(umr_mean)[-1,]
idh1 <- as.numeric(umr_mean_ma[,7]) 
idh2 <- as.numeric(umr_mean_ma[,8]) 
idh3 <- as.numeric(umr_mean_ma[,9]) 
no_per_number <- length(which((idh1>0.2)|(idh2>0.2)|(idh3>0.2))) #6505
no_per <- round(no_per_number/nrow(umr_mean),2) #0.26
no_no <- 1-no_per #0.74
no_po <- 1e-10
per_per <- 0.98
per_no <- 0.02
per_po <- 1e-10
po_per <- 1e-10
po_no <- 0.02
po_po <- 0.98
transMat <- rbind(c(per_per, per_no, per_po), c(no_per, no_no, no_po), c(po_per, po_no, po_po))


####RHmm get states of sites in reference UMR
#install.packages("nlme")
#install.packages("MASS")
#install.packages("RHmm")
library(nlme)
library(MASS)
library(RHmm)
###read methylation data from python-output files 
mutregion <- read.table("E:\\IDH\\R\\mutregion.txt")
wtregion <- read.table("E:\\IDH\\R\\wtregion.txt")

#dim(mutregion) #1777474       5
#dim(wtregion) #1763179       5
options(digits = 5)
mut_reg <- mutregion[,c(1,2,3,5)]
wt_reg <- wtregion[,c(1,2,3,5)]

###calculate states of sites
site_states <- function(gbm_reg){
        start = 1
        end = 1
        states = c()
        #Calculate the state of each ref_UMR separately
        for (i in 2:nrow(gbm_reg)){
          if((i!=nrow(gbm_reg))&(gbm_reg[i-1,1]==gbm_reg[i,1])&(gbm_reg[i-1,2]==gbm_reg[i,2])){
            end = end + 1
          }else{
            if(i==nrow(gbm_reg)){end = end +1}
            methy <- gbm_reg[start:end,4]
            dis <- distributionSet("NORMAL", means,sds)  #set normal distribution 
            initProb <- means
            hmm1 <- HMMSet(initProb, transMat, dis)  #set HMMSet
            VitPath <- viterbi(hmm1, methy)
            states <- c(states,VitPath$states)
            start = end +1
            end = end +1
          }
        }
        return(states)
}

mut_states <- site_states(mut_reg)
#mut_states[1288124:1288132] #2 1 1 1 1 1 1 1 1
wt_states <- site_states(wt_reg)
#wt_states[1288124:1288132] #2 2 2 2 2 2 2 2 2

###Replace isolated single and double sites with continuous ones, such as 121→111，1221→1111
con_states <- function(states1){
        ind <- c()
        for (i in 1:(length(states1)-3)){
          if(states1[i]==states1[i+3]){
            if((states1[i]!=states1[i+2])|(states1[i]!=states1[i+1])){ind<-c(ind,i)}
            states1[i+1]=states1[i]
            states1[i+2]=states1[i]
          }
        }
        return(states1)
}

mut_con_states <- con_states(mut_states)
wt_con_states <- con_states(wt_states)

###Put the States back on the CG site
mut_re <- mutregion
mut_re[,5] <- mut_con_states
wt_re <- wtregion
wt_re[,5] <- wt_con_states #chr,start,end,site,state

##output result to files
write.table(mut_re, "E:\\IDH\\R\\mutstates.txt",quote = FALSE,row.names = FALSE,sep = "\t",append = FALSE,col.names=FALSE)
write.table(wt_re, "E:\\IDH\\R\\wtstates.txt",quote = FALSE,row.names = FALSE,sep = "\t",append = FALSE,col.names=FALSE)

