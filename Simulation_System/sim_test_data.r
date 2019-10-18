library(MASS)
library(mvtnorm)
library(cluster)

# model1;

GenerateSpherePoints <- function(nrPoints,nrDim,center=rep(0,nrDim),r=1){
  #generate the polar coordinates!
  x <-  matrix(runif(nrPoints*nrDim,-pi,pi),ncol=nrDim)
  x[,nrDim] <- x[,nrDim]/2
  #recalculate them to cartesians
  sin.x <- sin(x)
  cos.x <- cos(x)
  cos.x[,nrDim] <- 1  # see the formula for n.spheres
  
  y <- sapply(1:nrDim, function(i){
    if(i==1){
      cos.x[,1]
    } else {
      cos.x[,i]*apply(sin.x[,1:(i-1),drop=F],1,prod)
    }
  })*sqrt(runif(nrPoints,0,r^2))
  
  y <-  as.data.frame(
    t(apply(y,1,'+',center))
  )
  
  names(y) <- make.names(seq_len(nrDim))
  y
}


x1<-GenerateSpherePoints(200,10)
plot(x1)
pairs(x1)

#loop;
x1.ch<-NULL
x1.kl<-NULL
x1.hart<-NULL
x1.sil<-NULL
x1.gap<-NULL
ind<-NULL
for (i in 1:50) {
  set.seed(i+50)
  x1<-GenerateSpherePoints(200,10)
  x1.ch[i]<-ch(x1)
  x1.kl[i]<-kl(x1)
  x1.hart[i]<-hart(x1)
  x1.sil[i]<-sil(x1)
  x1.gap[i]<-gap(x1)
  in1<-rbind(x1.gap,x1.sil,x1.ch,x1.kl,x1.hart)
  ind1<-apply(ind1,1,table)
}
print(ind1)

#####################################################################

#model2;

mu1<-c(0,0)
mu2<-c(0,5)
mu3<-c(5,-3)
n1=25
n2=25
n3=50
sigma=matrix(c(1,0.5,0.5,1),2,2)
x2a=mvrnorm(n1, mu = mu1, Sigma = sigma)
x2b=mvrnorm(n2, mu = mu2, Sigma = sigma)
x2c=mvrnorm(n3, mu = mu3, Sigma = sigma)
x2<-rbind(x2a,x2b,x2c)
plot(x2)


#loop;

x2a=matrix()
x2b=matrix()
x2c=matrix()
x2=matrix()
mu1<-c(0,0)
mu2<-c(0,5)
mu3<-c(5,-3)
n1=25
n2=25
n3=50
sigma=matrix(c(1,0.5,0.5,1),2,2)
x2.ch<-NULL
x2.kl<-NULL
x2.hart<-NULL
x2.sil<-NULL
x2.gap<-NULL
ind2<-NULL
for (i in 1:50) {
  set.seed(i+100)
  x2a=mvrnorm(n1, mu = mu1, Sigma = sigma)
  x2b=mvrnorm(n2, mu = mu2, Sigma = sigma)
  x2c=mvrnorm(n3, mu = mu3, Sigma = sigma)
  x2<-rbind(x2a,x2b,x2c)
  
  x2.ch[i]<-ch(x2)
  x2.kl[i]<-kl(x2)
  x2.hart[i]<-hart(x2)
  x2.sil[i]<-sil(x2)
  x2.gap[i]<-gap(x2)
  in2<-data.frame(x2.gap,x2.sil,x2.ch,x2.kl,x2.hart)
  ind2<-apply(in2,2,table)
}
print(ind2)

#####################################################################

#model3;
set.seed(123)
(n1=ifelse(rbinom(1,1,0.5)==0,25,50))
(n2=ifelse(rbinom(1,1,0.5)==0,25,50))
(n3=ifelse(rbinom(1,1,0.5)==0,25,50))
(n4=ifelse(rbinom(1,1,0.5)==0,25,50))
n=n1+n2+n3+n4
X=matrix(0,n,10)
colnames(X)=paste("X",1:ncol(X),sep = "")
p=3
set.seed(123)
mu11=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
mu1=cbind(t(mu11),t(rep(0,7)))
mu22=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
mu2=cbind(t(mu22),t(rep(0,7)))
mu33=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
mu3=cbind(t(mu33),t(rep(0,7)))
mu44=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
mu4=cbind(t(mu44),t(rep(0,7)))
x3a=rmvnorm(n1,mu1,diag(10))
x3b=rmvnorm(n2,mu2,diag(10))
x3c=rmvnorm(n3,mu3,diag(10))
x3d=rmvnorm(n4,mu4,diag(10))
x3=rbind(x3a,x3b,x3c,x3d)
pairs(x3)

#loop;

n1<-NULL
n2<-NULL
n3<-NULL
n4<-NULL
n<-NULL
mu11<-NULL
mu1<-NULL
mu22<-NULL
mu2<-NULL
mu33<-NULL
mu3<-NULL
mu44<-NULL
mu4<-NULL
x3a=matrix()
x3b=matrix()
x3c=matrix()
x3d=matrix()
x3=matrix()
x3.ch<-NULL
x3.kl<-NULL
x3.hart<-NULL
x3.sil<-NULL
x3.gap<-NULL
in3<-NULL
ind3<-NULL
for (i in 1:50) {
  set.seed(i+50)
  n1=ifelse(rbinom(1,1,0.5)==0,25,50)
  n2=ifelse(rbinom(1,1,0.5)==0,25,50)
  n3=ifelse(rbinom(1,1,0.5)==0,25,50)
  n4=ifelse(rbinom(1,1,0.5)==0,25,50)
  n=n1+n2+n3+n4
  p=3
  
  mu11=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
  mu1=cbind(t(mu11),t(rep(0,7)))
  mu22=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
  mu2=cbind(t(mu22),t(rep(0,7)))
  mu33=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
  mu3=cbind(t(mu33),t(rep(0,7)))
  mu44=sample(rmvnorm(3,rep(0,p),25*diag(p)),3)
  mu4=cbind(t(mu44),t(rep(0,7)))
  x3a=rmvnorm(n1,mu1,diag(10))
  x3b=rmvnorm(n2,mu2,diag(10))
  x3c=rmvnorm(n3,mu3,diag(10))
  x3d=rmvnorm(n4,mu4,diag(10))
  x3=rbind(x3a,x3b,x3c,x3d)
  
  x3.ch[i]<-ch(x3)
  x3.kl[i]<-kl(x3)
  x3.hart[i]<-hart(x3)
  x3.sil[i]<-sil(x3)
  x3.gap[i]<-gap(x3)
  in3<-data.frame(x3.gap,x3.sil,x3.ch,x3.kl,x3.hart)
  ind3<-apply(in3,2,table)
}
print(ind3)


##################################################################
#model4;
set.seed(1234)
(n1=ifelse(rbinom(1,1,0.5)==0,25,50))
(n2=ifelse(rbinom(1,1,0.5)==0,25,50))
(n3=ifelse(rbinom(1,1,0.5)==0,25,50))
(n4=ifelse(rbinom(1,1,0.5)==0,25,50))
n=n1+n2+n3+n4
X=matrix(0,n,10)
colnames(X)=paste("X",1:ncol(X),sep = "")
p=10
set.seed(123)
mu1=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
mu2=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
mu3=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
mu4=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
x4a=rmvnorm(n1,mu1,diag(10))
x4b=rmvnorm(n2,mu2,diag(10))
x4c=rmvnorm(n3,mu3,diag(10))
x4d=rmvnorm(n4,mu4,diag(10))
x4=rbind(x4a,x4b,x4c,x4d)
pairs(x4)

#loop;
n1<-NULL
n2<-NULL
n3<-NULL
n4<-NULL
n<-NULL
mu1<-NULL
mu2<-NULL
mu3<-NULL
mu4<-NULL
x4a=matrix()
x4b=matrix()
x4c=matrix()
x4d=matrix()
x4=matrix()
x4.ch<-NULL
x4.kl<-NULL
x4.hart<-NULL
x4.sil<-NULL
x4.gap<-NULL
x4.gapPC<-NULL
in4<-NULL
ind4<-NULL
for (i in 1:50) {
  set.seed(i+50)
  n1=ifelse(rbinom(1,1,0.5)==0,25,50)
  n2=ifelse(rbinom(1,1,0.5)==0,25,50)
  n3=ifelse(rbinom(1,1,0.5)==0,25,50)
  n4=ifelse(rbinom(1,1,0.5)==0,25,50)
  n=n1+n2+n3+n4
  p=10
  
  mu1=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
  mu2=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
  mu3=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
  mu4=sample(rmvnorm(3,rep(0,p),3.6*diag(p)),10)
  x4a=rmvnorm(n1,mu1,diag(10))
  x4b=rmvnorm(n2,mu2,diag(10))
  x4c=rmvnorm(n3,mu3,diag(10))
  x4d=rmvnorm(n4,mu4,diag(10))
  x4=rbind(x4a,x4b,x4c,x4d)
  
 
  x4.ch[i]<-ch(x4)
  x4.kl[i]<-kl(x4)
  x4.hart[i]<-hart(x4)
  x4.sil[i]<-sil(x4)
  x4.gap[i]<-GAP(x4)
  x4.gapPC[i]<-gapPC(x4)
  
  in4<-data.frame(x4.gap,x4.gapPC,x4.sil,x4.ch,x4.kl,x4.hart)
  ind4<-apply(in4,2,table)
}
print(ind4)

#################################################################
#model5;
n=100

