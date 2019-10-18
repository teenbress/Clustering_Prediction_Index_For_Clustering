setwd("C:/Users/Qiao Yu/Desktop/data")
library(iCluster)
data(breast.chr17)
summary(breast.chr17)
attach(breast.chr17)
dim(mRNA.data)
dim(DNA.data)

##separate clustering by heatmap;
library(gplots)

my_palette <- colorRampPalette(c("darkblue", "white","red"))(n =300)

heatmap.2(mRNA.data,col = my_palette,scale="none",breaks = seq(-1.5,1.5,0.01),trace = "none",dendrogram="row",margins=c(4,8))

my_palette <- colorRampPalette(c("darkblue", "white","red"))(n =180)

heatmap.2(DNA.data,col = my_palette,scale="none",breaks = seq(-0.9,0.9,0.01),trace = "none",dendrogram="row",margins=c(4,8))

## fig.2 B
fit=iCluster(breast.chr17, k=4, lambda=c(0.2,0.2))
plotiCluster(fit=fit, label=rownames(breast.chr17[[2]]))
compute.pod(fit)
order.idx <- order(fit$cluster)
mRNA.order <- breast.chr17$mRNA.data[order.idx,]
DNA.order<-DNA.data[order.idx,]

my_palette <- colorRampPalette(c("darkblue", "white","red"))(n =300)

heatmap.2(mRNA.order,col = my_palette,scale="none",breaks = seq(-1.5,1.5,0.01),trace = "none",dendrogram="none",Rowv=FALSE,Colv=FALSE,margins = c(4,8))

my_palette <- colorRampPalette(c("darkblue", "white","red"))(n =180)

heatmap.2(DNA.order,col = my_palette,scale="none",breaks = seq(-0.9,0.9,0.01),trace = "none",dendrogram="none",Rowv=FALSE,Colv=FALSE,margins=c(4,8))



# k-means clustering
cl <- kmeans(mRNA.data, centers=4)

# use order to sort a vector named $cluster in list cl
# it returns an index
(order.ikx <- order(cl$cluster))

# use the index to sort the original matrix by row
# note: row 1: [1,] / column 1: [,1]
df.m.order <- breast.chr17$mRNA.data[order.ikx,]
heatmap(df.m.order)


##the real data: the mythelation part d1;the gene exp part d2;
a<- read.table("S5.txt", header=T, sep="\t")
aa<-a[,-1]
sdd<-apply(aa,1,sd)
sde<-sapply(a[1:10,],sd)
length(which(sdd<=0.313))#exclude 1164 CpG sites
methy<-a[which(sdd>0.313),]
b<- read.table("S1.txt", header=T, sep="\t")
d2<-b[,-1]
d1<-methy
dim(d1)
dim(d2)

asvd<-svd(scale(d2[,-1]))
UE <- asvd$u %*% diag(asvd$d)
PC1 <- UE[,1]
PC2 <- UE[,2]
variance <- asvd$d^2 / sum(asvd$d^2)
v1 <- paste0("variance: ",signif(variance[1] * 100,3), "%")
v2 <- paste0("variance: ",signif(variance[2] * 100,3), "%")
plot(PC1, PC2, col=as.numeric(iris$Species),pch=19, xlab=v1, ylab=v2)