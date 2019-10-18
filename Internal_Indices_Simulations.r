# Internal Indices fuctions and clustering simulation
# Including£º
	
#   Calinski and Harabasz (ch)	
#   Krzanowski and Lai Index (kl)
#   The silhouette width method (sil)
#   Gap
#   Gap (pc)	
#   Hartigan(hart)	


library(factoextra)
library(cluster)
library(NbClust)
library(cclust)
library(mvtnorm)

# gap, sil,ch, kl, hart for PAM clustering; 
rm(ch)
rm(kl)
rm(sil)
rm(hart)

ch <- function(x)NbClust(x, distance = "euclidean", min.nc = 2,max.nc = 10, method = "ward.D", index ="ch")$Best.nc[1]

kl <- function(x)NbClust(x, distance = "euclidean", min.nc = 2,max.nc = 10, method = "complete", index ="kl")$Best.nc[1]

sil <- function(x)NbClust(x, distance = "euclidean", min.nc = 2,max.nc = 10, method = "complete", index ="sil")$Best.nc[1]

GAP<- function(x){
 gapn<-clusGap(x, spaceH0 = "original",FUN = pam, K.max = 10, B = 20,verbose = FALSE)
 gap.ind<- maxSE(gapn$Tab[, "gap"], gapn$Tab[, "SE.sim"], method="Tibs2001SEmax")
 return(gap.ind)
 }

gapPC<- function(x){
  gappc<-clusGap(x, spaceH0 = "scaledPCA",FUN = pam, K.max = 10, B = 20,verbose = FALSE)
  gapPC.ind<- maxSE(gappc$Tab[, "gap"], gappc$Tab[, "SE.sim"], method="Tibs2001SEmax")
  return(gapPC.ind)
}

hart<- function(x, k){
  n <- nrow(x)
  clall<-function(hc,nc)
  {
    cl1 <- pam(hc, k=nc)$clustering
    cl2 <- pam(hc, k=nc+1)$clustering
    clall <- cbind(cl1, cl2)
    return(clall)
  }
  
  gss <- function(x, cl) 
  {
    n <- length(cl)
    k <- max(cl)
    centers<-pam(x,k)$medoids
    withins <- rep(0, k)
    x.2 <- (x - centers[cl, ])^2
    for (i in 1:k) {
      withins[i] <- sum(x.2[cl == i, ])
    }
    wgss <- sum(withins)
    results <- list(wgss=wgss, centers=centers)
    return(results)
  }
  
  
  HART<-NULL
  for(i in 1:k){
    call<-clall(x,i)
    HART[i] <- (gss(x, call[, 1])$wgss/gss(x, call[, 2])$wgss - 1.039) * (n - k - 1)
  }
  num.hart<-min(which(HART<=10))
  return(num.hart)
}

#hart(c,k=10)

#loop for t= t times;

rm(indexfun)
indexfun <- function(x,...){ 
  c( gap = GAP(x), gapPC = gapPC(x), sil = sil(x), ch = ch(x), kl = kl(x), hart = hart(x,k=10)) 
}

rm(result)
result<-function(t){
  a = matrix(NA,nrow=t,ncol=6)
  
  colnames(a) <- c("gap","gapPC","sil","ch","kl","hart")
  
  simu1 = simu2 = simu3 = simu4 = simu5 = simu6 = simu7 = simu8 = a
  
  for (i in 1:t){
    d1 <- sim_a()
    d2 <- sim_b()$data
    d3 <- sim_c()$data
    d4 <- sim_d()$data
    d5 <- sim_e()$data
    d6 <- sim_f()$data
    d7 <- sim_g()$data
    d8 <- sim_h()$data
    
    simu1[i,] <- indexfun(d1)
    simu2[i,] <- indexfun(d2)
    simu3[i,] <- indexfun(d3)
    simu4[i,] <- indexfun(d4)
    simu5[i,] <- indexfun(d5)
    simu6[i,] <- indexfun(d6)
    simu7[i,] <- indexfun(d7)
    simu8[i,] <- indexfun(d8)
    simu <- list(simu1, simu2, simu3, simu4, simu5, simu6, simu7, simu8)
  }
  print(simu)
}


a <- result(t = 100)
(final <- lapply(a, function(x)apply(x,2,table)))
write.csv(a,"a.csv")




