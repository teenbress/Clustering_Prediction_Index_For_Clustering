

library(sfsmisc)
library(robustbase)

##########################################################################
# External Index;

rm(EI)
EI <- function(x, k, B=20, n.tr=floor(2/3*nrow(x)),
               clust= function(x,k) pam(x,k, cluster.only=TRUE),
               class= diagDA,
               mit = median, ext.ind="FM", kor=FALSE, part=FALSE)
{n <- nrow(x)
s <- numeric(B)
n.test <- n-n.tr
index <- sapply(1:B,function(i) sample(n,n.tr))

for(b in 1:B){
  train <- x[index[,b],]
  test <- x[-index[,b],]
  test <- na.exclude(test)
  klassen.train <- clust(x=train,k=k)
  test.pred <- class(ls=train,cll=klassen.train,ts=test)
  klassen.test <- clust(x=test,k=k)
  tab <- table(test.pred,klassen.test)
  s.row <- rowSums(tab)
  s.col <- colSums(tab)
  if(kor){
    if(ext.ind=="FM"){
      ra <- sum(tab*(tab-1)/2)
      prod <- sum(s.row*(s.row-1)/2)*sum(s.col*(s.col-1)/2)
      E <- prod/(n.test*(n.test-1)/2)
      s[b] <- (ra-E)/(sqrt(prod)-E)
    }
    if(ext.ind=="Rand"){
      ra <- sum(tab*(tab-1)/2)
      sum.row <- sum(s.row*(s.row-1)/2)
      sum.col <- sum(s.col*(s.col-1)/2)
      E <- sum.row*sum.col/(n.test*(n.test-1)/2)
      s[b] <- (ra - E)/((sum.row+sum.col)/2-E)
    }
  }else{
    if(ext.ind=="FM"){
      Z <- sum(tab^2)
      s[b] <- 1/2*(Z-nrow(test))/(sum(s.row*(s.row-1)/2)*
                                    sum(s.col*(s.col-1)/2))^(1/2)
    }
    if(ext.ind=="Rand"){
      ra <- sum(tab*(tab-1)/2)
      n.2 <- n.test*(n.test-1)/2
      sum.row <- sum(s.row*(s.row-1)/2)
      sum.col <- sum(s.col*(s.col-1)/2)
      s[b] <- 1 + (2*ra - (sum.row+sum.col))/n.2
    }
  }
}
res <- mit(s)
list(Index=res, S=s, part = index)
}  

# End of the function EI;
# Classifier DLDA;

#########################################################################
# Function to estimate the number of clusters;

rm(k.opt)

k.opt <- function(p.k, d.k, p.max, d.min, name=FALSE)
{
  ausw <- function(j,i)
  {
    d.k <- d.k.org
    low <- (j-1)*n.d.k + 1
    up <- j*n.d.k
    if(any(K.min[low:up,i])){
      d.k[!K.min[low:up,i]] <- 0
      which.max(d.k) + 1
    }else 1
  }
  d.k.org <- d.k
  n.d.k <- length(d.k)
  K.min <- sapply(p.max, function(p.max){
    sapply(d.min, function(d.min) p.k <= p.max & d.k >= d.min)})
  if(!is.matrix(K.min)){
    res <- if(any(K.min)) 2 else 1
  }else{
    res <- sapply(1:length(p.max), function(i){
      sapply(1:length(d.min), function(j) ausw(j,i))
    })
  }
  if(name){
    dimnames(res) <- list(as.character(d.min),as.character(p.max))
  }
  res ##  represent the d.min, columns the d.max values;
}

unif.null <- function(x, rx, mom, n.set, d,n)
{
  matrix(runif(n.set,rx[1,],rx[2,]),ncol=d, byrow=TRUE)
}

#############################################################################
rm(clest)

clest <- function(x, M, p.max=0.05, d.min=0.05, B=20, n0=25,
                  n.tr=floor(2/3*nrow(x)),
                  clust = function(x, k) pam(x=x, k=k, cluster.only=TRUE),
                  class = diagDA,
                  null.model = unif.null, ext.ind = "FM",
                  mit = function(x) huberM(x=x)$mu, kor=FALSE,a=.05)
{
  if(!is.matrix(x)){
    x <- data.matrix(x)
  }
  n <- nrow(x)
  d <- ncol(x)
  n.na <- n-nrow(na.exclude(x))
  name <- as.character(seq(2,M))
  
  t <- numeric(M-1)
  S <- matrix(0,B,M-1)
  t.ref <- matrix(0,nrow=M-1,ncol=n0)
  ##Max. And Min. For the Unif null model;
  rx <- apply(x,2,range,na.rm=TRUE)
  
  df <- apply(x,2,ecdf)
   n.set <- d*n
   S.ref <- array(0,dim=c(B,M-1,n0),dimnames=list(NULL,name,NULL))
  ## save null sets;
  null.sets <- array(0,dim=c(n,d,n0))
  for(j in 1:n0){
    null.sets[,,j] <- null.model(x=x,rx=rx,n.set=n.set, d=d)
    {matrix(runif(n.set,rx[1,],rx[2,]),ncol=d, byrow=TRUE)}
  }
  stat <- EI(x=x,k=2,B=B,n.tr=n.tr,clust=clust,class=class,mit=mit,
             ext.ind=ext.ind, kor=kor)
  t[1] <- stat$Index
  S[,1] <- stat$S
  part <- stat$part
  part0 <- array(0,dim=c(n.tr,B,n0))
  for(i in 1:n0){
    res.null <- EI(x=null.sets[,,i],k=2,B=B,n.tr=n.tr,clust=clust,
                   class=class, mit=mit, ext.ind=ext.ind, kor=kor)
    S.ref[,1,i] <- res.null$S
    t.ref[1,i] <- res.null$Index
    part0[,,i] <- res.null$part
  }
  for(k in 3:M){
    stat <- EI(x=x,k=k,B=B,n.tr=n.tr,clust=clust,class=class,mit=mit,
               ext.ind=ext.ind, kor=kor,part=part)
    t[k-1] <- stat$Index
    S[,k-1] <- stat$S
    for(i in 1:n0){
      res.null <- EI(x=null.sets[,,i],k=k,B=B,n.tr=n.tr,clust=clust,
                     class=class, mit=mit, ext.ind=ext.ind, kor=kor,
                     part=part0[,,i])
      S.ref[,k-1,i] <- res.null$S
      t.ref[k-1,i] <- res.null$Index
    }
  }
  t.ref.0 <- rowMeans(t.ref)
  p.k <- rowMeans(t.ref>=t)
  d.k <- t-t.ref.0
  K.t <- which(t == max(t)) + 1
  K <- k.opt(p.k=p.k,d.k=d.k,p.max=p.max,d.min=d.min)
  temp <- d.k
  temp[K-1] <- 0
  K2 <- k.opt(p.k=p.k,d.k=temp,p.max=p.max,d.min=d.min)
  names(t) <- names(t.ref.0) <- names(d.k) <- names(p.k) <- name
  dimnames(S) <- list(NULL,name)
  dimnames(t.ref) <- list(name,NULL)
  Liste <- list(call= match.call(), K = K, K2=K2, K.t = K.t, p.k = p.k,
                p.max = p.max, d.k = d.k, n.tr=n.tr,
                t=t, S=S, t.ref = t.ref, t.ref.0 = t.ref.0, S.ref=S.ref)
  ## Liste$KIs <- confint(x=Liste)
  class(Liste) <- "clest"
  Liste
}

#clest test;
# clest(d,M=6)$K
rm(clust.result)
# clest loops for t times ########
clust.result<-function(t){
  
  a=matrix(NA,nrow=t,ncol=1)
  
  colnames(a)<-c("Clest")
  
  simu1=simu2=simu3=simu4=simu5=simu6=simu7=simu8=a
  
  for (i in 1:t){
    d1<-sim_a()
    d2<-sim_b()$data
    d3<-sim_c()$data
    d4<-sim_d()$data
    d5<-sim_e()$data
    d6<-sim_f()$data
    d7<-sim_g()$data
    d8<-sim_h()$data
    simu1[i,]<- clest(d1,M=6)$K
    simu2[i,]<- clest(d2,M=6)$K
    simu3[i,]<- clest(d3,M=6)$K
    simu4[i,]<- clest(d4,M=6)$K
    simu5[i,]<- clest(d5,M=6)$K
    simu6[i,]<- clest(d6,M=6)$K
    simu7[i,]<- clest(d7,M=6)$K
    simu8[i,]<- clest(d8,M=6)$K
  }
  simu<-cbind(simu1,simu2,simu3,simu4,simu5,simu6,simu7,simu8)
  res<-apply(simu,2,table)
  return(simu)
  print(res)
}

clust.result(t=100)

a<-read.table("a.txt",sep = "")
a<-t(a[,-1])
apply(a,1,table)

      