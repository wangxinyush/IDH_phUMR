cgi_methy_file <- read.table("CGI_brain_methy_matrix.txt", header=T, sep="\t", fill=T,
                  check.names=F, quote="\"", dec = ".", comment.char = "")
dim(cgi_methy_file) #27718 92
cgi_methy_file<-na.omit(cgi_methy_file)
dim(cgi_methy_file)  #6848 92
cgi_methy_file[1:10,1:10]

pca2 <- prcomp(cgi_methy_file[,4:92], scale. =F, center=F )

Tissue<-colnames(cgi_methy_file)[4:92]
NeuN<-colnames(cgi_methy_file)[4:92]

df<-data.frame(samples=colnames(cgi_methy_file)[4:92],
               Tissue=Tissue,
               NeuN=NeuN,
               pc1=pca2$rotation[,1], pc2=pca2$rotation[,2]
)
head(df)
df[,2:3]

colors<-c("#E41A1C","#FF7F00","#800101","#4DAF4A", "#377EB8","#984EA3")
ggplot(df, aes(x= pc1, y= pc2,shape=NeuN))+
  geom_point(aes(color =Tissue),size=2)+
  scale_colour_manual(values=alpha(c(colors,0.2)))+
  theme_bw()+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size = 20),
        axis.text=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))

ggplot(df, aes(x= pc1, y= pc2))+
  geom_point(aes(color =NeuN),size=2)+
  theme_classic()+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size = 20),
        axis.text=element_text(size=20),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))

#pca plot
pca1 <- prcomp(t(cgi_methy_file[,4:92]), scale. =TRUE, center=TRUE)

screeplot(pca1,type='line',npcs = length(pca1$sdev), lwd=2, col=blue)
plot(pca1$x[,1], pca1$x[,2], col="blue", pch=19, xlab="PC1", ylab="PC2", cex.lab=1.3)
points(pca1$x[c(8, 9, 10),1], pca1$x[c(8, 9, 10),2], col="red", pch=19)
#text(pca1$x[,1], pca1$x[,2], cex=1, labels = c("A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2", "F2", "G2", "H2") )
text(pca1$x[,1], pca1$x[,2], cex=1, labels = c("A1", "B1", "C1", "D1", "E1", "G1", "H1", "B2", "E2", "G2", "H2") )



