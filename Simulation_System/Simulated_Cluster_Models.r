####################################################
#The simulation section datasets
####################################################
library(mvtnorm)

#Model 1

sim_a <- function(nclus=1, ncov=10, clustersize=200,
                              tmin=0, tmax=1, sepfactor=10)
{
  simdat <- matrix(0, nrow=nclus*clustersize, ncol=ncov)
  offset <- 1
  for (i in 1:nclus) {
    for (j in 1:ncov) {
      simdat[offset:(offset + clustersize - 1), j] <- 
        runif(clustersize,
              min=((i-1)*sepfactor + tmin),
              max=((i-1)*sepfactor + tmax))
    }
    offset <- offset + clustersize
  }
  return(simdat)
}
####################################################
#Model 2

sim_b <- function(
  counts = c(25, 50, 25),
  centers = matrix(c(0, 0, 0, 5, 5, -3), ncol = 2, byrow = T),
  mu = 0, sdev = 1)
{
  if (length(counts) != nrow(centers))
    stop("Counts must have same length as number of rows.")
  
  eps <- rnorm(sum(counts) * ncol(centers), mean = mu, sd = sdev)
  
  mat <- Map(function(row_i, cnt)
  {
    cur_row <- centers[row_i,]
    matrix(rep(cur_row, cnt), ncol = ncol(centers), byrow = T)
  }, 1:nrow(centers), counts)
  
  mat <- do.call(rbind, mat)
  mat <- mat[sample.int(nrow(mat)),] + eps
  
  list(data = mat, nclust = nrow(centers), ndim = ncol(centers))
}
########################################################################

#Model 3
#' Replicate centers

rep_centers <- function(centers, counts)
{
  mat <- Map(function(row_i, cnt)
  {
    cur_row <- centers[row_i,]
    matrix(rep(cur_row, cnt), ncol = ncol(centers), byrow = T)
  }, 1:nrow(centers), counts)
  
  do.call(rbind, mat)
}

#' Generate random multivariate normal centers

gen_mvnorm_centers <- function(nclust, mu, sigma, min_dist, max_iter = 1000)
{
  stopifnot(length(mu) == nrow(sigma) && nrow(sigma) == ncol(sigma))
  
  cur_iter <- 0
  while(cur_iter < max_iter) {
    centers <- MASS::mvrnorm(nclust, mu, sigma)
    
    if (all(dist(centers) >= min_dist))
      return(centers)
    
    cur_iter <- cur_iter + 1
  }
  
  stop("Couldn't generate centers with min_dist by max_iter.")
}

sim_c <- function(
  nclust = 4,
  nobs = c(25, 50),
  mu = rep(0, 4),
  sigma = diag(5, 4)
)
{
  clust_size <- sample(nobs, nclust, replace = T)
  centers <- gen_mvnorm_centers(nclust, mu, sigma, 1)
  sim_vals <- rep_centers(centers, clust_size)
  eps <- rnorm(length(sim_vals))
  sim_vals <- sim_vals + eps
  sim_vals <- sim_vals[sample.int(nrow(sim_vals)),]
  
  list(data = sim_vals, clust_size = clust_size, nclust = length(clust_size),
       ndim = ncol(sim_vals))
}

#############################################################
#model 4
sim_d <- function(nobs = c(25, 50))
{
  sim_c(nclust = 4,
        nobs = nobs,
        mu = rep(0, 10),
        sigma = diag(3.6, 10))
}

##########################################################

#model5

sim_e <- function(nsamp=100, rng = c(-0.5, 0.5))
{
  clust_1 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
                    ncol = 3)
  clust_2 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
                    ncol = 3) + 10
  
  clust_1 <- clust_1 + rnorm(length(clust_1), mean = 0, sd = 0.1)
  clust_2 <- clust_2 + rnorm(length(clust_2), mean = 0, sd = 0.1)
  
  res <- rbind(clust_1, clust_2)
  res <- res[sample.int(nrow(res)),]
  
  list(data = res, nclust = 2, ndim = ncol(clust_1))
}


###############################################################
#model 6

sim_f <- function(nsamp=100, rng = c(-0.5, 0.5))
{
  clust_1 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
                    ncol = 3)
  clust_2 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
                    ncol = 3) + 10
  
  clust_1 <- clust_1 + rnorm(length(clust_1), mean = 0, sd = 0.1)
  clust_2 <- clust_2 + rnorm(length(clust_2), mean = 0, sd = 0.1)
  res1 <- rbind(clust_1, clust_2)
  
  v=c(4:10)
  noise<-function(x)rnorm(nrow(res1), mean = 0, sd = x)
  res2<-sapply(v,noise)
  
  res<-cbind(res1,res2)
  res <- res[sample.int(nrow(res)),]
  
  list(data = res, nclust = 2, ndim = ncol(clust_1))
}


################################################################
#model 7

sim_g  <- function(nsamp=50)
{
  clust_1 <- matrix(rnorm(nsamp, mean = 0, sd = 1),ncol=1)
  clust_2 <- matrix(rnorm(nsamp, mean = 2.5, sd = 1),ncol=1)
  res1 <- rbind(clust_1, clust_2)
  res2 <- rmvnorm(nsamp*2,mean=rep(0,9),sigma=diag(9))
 
  res <- cbind(res1,res2)
  res <- res[sample.int(nrow(res)),]
  
  list(data =res , nclust = 2, ndim = ncol(clust_1))
}

################################################################
#model 8
sim_h  <- function(nsamp=50)
{ 
  sig=matrix(rep(0.5,9),ncol=3)
  diag(sig)=1
  clust_1 <- rmvnorm(nsamp,mean=c(0,0,0),sigma=sig)
  clust_2 <- rmvnorm(nsamp,mean=c(2,-2,2),sigma=sig)
  clust_3 <- rmvnorm(nsamp,mean=c(-2,2,-2),sigma=sig)
  res1 <- rbind(clust_1, clust_2,clust_3)
  res2 <- rmvnorm(nsamp*3,mean=rep(0,10),sigma=diag(10))
  
  res <- cbind(res1,res2)
  res <- res[sample.int(nrow(res)),]
  
  list(data =res , nclust = 2, ndim = ncol(clust_1))
}






#############################################################

#Simulate 8 model dataset one time 

multi.fun <- function(...) {
     list( sim_a(), sim_b()$data,  sim_c()$data, sim_d()$data, sim_e()$data,sim_f()$data,
           sim_g()$data, sim_h()$data)
}
d<-multi.fun()


######################################################################################
#loop for t= t times,moved to simulate index.r;


indexfun<-function(x,...){ 
  c(gap=GAP(x),gapPC=gapPC(x),sil=sil(x),ch=ch(x),kl=kl(x),hart=hart(x)) 
}



result<-function(t){
  
  
  a=matrix(NA,nrow=t,ncol=6)
  
  colnames(a)<-c("gap","gapPC","sil","ch","kl","hart")
  
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
    simu1[i,]<- indexfun(d1)
    simu2[i,]<- indexfun(d2)
    simu3[i,]<- indexfun(d3)
    simu4[i,]<- indexfun(d4)
    simu5[i,]<- indexfun(d5)
    simu6[i,]<- indexfun(d6)
    simu7[i,]<- indexfun(d7)
    simu8[i,]<- indexfun(d8)
    simu<-list(simu1,simu2,simu3,simu4,simu5,simu6,simu7,simu8)
  }
  print(simu)
}
a<-result(t=50)
(final<-lapply(a, function(x)apply(x,2,table)))





