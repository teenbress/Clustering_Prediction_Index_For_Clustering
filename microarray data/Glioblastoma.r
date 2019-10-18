library(iClusterPlus)
library(GenomicRanges)
library(gplots)
library(lattice)
data(gbm)
dim(gbm.mut)
dim(gbm.exp)
dim(gbm.seg)
dim(variation.hg18.v10.nov.2010)
gbm.cn=CNregions(seg=gbm.seg,epsilon=0,adaptive=FALSE,rmCNV=TRUE,cnv=variation.hg18.v10.nov.2010[,3:5],frac.overlap=0.5, rmSmallseg=TRUE,nProbes=5)
dim(gbm.cn)
mut.rate=apply(gbm.mut,2,mean)
gbm.mut2 = gbm.mut[,which(mut.rate>0.02)]
dim(gbm.mut2)





## clustering
n1 = 20
n2 = 20
n3 = 20
n = n1+n2+n3
p = 5
q = 100

x = NULL
x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[1]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[2]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[3]] = cbind(xa,xb)


x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[4]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[5]] = cbind(xa,xb)

method = c('lasso', 'enet', 'flasso', 'glasso', 'gflasso')  
lambda=alist()
lambda[[1]] = 30
lambda[[2]] = c(20,1)
lambda[[3]] = c(20,20)
lambda[[4]] = 30
lambda[[5]] = c(30,20)

chr=c(rep(1,10),rep(2,(p+q)-10)) 
fit = iCluster(x, k=3, lambda)
fit2 = iCluster2(x, K=3, lambda, method=method, chr=chr, maxiter=20,epsilon=1e-8)


par(mfrow=c(1,1),mar=c(4,4,1,1))
for(i in 1:5){
  barplot(fit2$beta[[i]][,1])
}

