rm(list = ls()) 
setwd("E:/ECNU/Reaserch/matrix2sample_test/code_new/code")
source("function_set.r")
library(MASS) 
library(mvtnorm)
library(energy)    ## energy test
library(maotai)    ##Maximum Mean Discrepanc test
library(MixMatrix) ##generate matrix-variate normal/t distribution
library(lava)      ## trace operator
library(flipr)
library(gTests)
library(Ball)
library(spdep)
library(highmean)
library(ade4)
###### ##### ##### ##### ##### ##### 
###### ##### Simulation Setting ##### 
###### ##### ##### ##### ##### ##### ##### 
P <- 25 # number of row
Q <- 15    # number of column
rho<- 0.1
U <-  generate_sig(P,rho)
V <-  generate_sig(Q,rho)



###### ##### ##### ##### ##### ##### 

example=21
max_iter=8
N1=N2=20   
B <- 399    #permute times
nsim<- 500       #number of repetition
N_seq=rep(0, max_iter)
power=data.frame(Detc=numeric(), MMD=numeric(), Energy=numeric(),
                 Graph=numeric(),RGtest=numeric(),Ball=numeric(),
                 CQ=numeric(),Kim=numeric())


for (i in 1:1){
  
  
  rho2=rho
  U2 <-  generate_sig(P,rho2)
  V2 <-  generate_sig(Q,rho2)
  pvalues<-data.frame(Detc=numeric(), MMD=numeric(), Energy=numeric(),
                      Graph=numeric(),RGtest=numeric(),Ball=numeric(),
                      CQ=numeric(),Kim=numeric())
  
  for(loop in 1:nsim){
    
    ###### ##### ##### ##### ##### ##### 
    x <- rep(0,P*Q*N1)
    dim(x) <- c(P,Q,N1)
    y <- x
    
    if (example==21){
      for (j in 1:N1){
        x[,,j] <- rmatrixnorm(n = 1, mean=matrix(0,P,Q),U=U,V=V)
        y[,,j] <- rmatrixnorm(n = 1, mean=matrix(0,P,Q),U=U,V=V)
        #y[,,j] <- rmatrixnorm(n = 1, mean=matrix(0,P,Q),U=U,V=V2)
        #y[,,j] <- rmatrixnorm(n = 1, mean=matrix(0,P,Q),U=U2,V=V2)
      }}
    
    
    if (example==22){
      for (j in 1:N1){
        df=1
        x[,,j] <- rmatrixt(n=1,df,mean=matrix(0,P,Q), U, V)
        y[,,j] <-rmatrixt(n=1,df,mean=matrix(0,P,Q), U2, V)
        #y[,,j] <- rmatrixt(n=1,df,mean=matrix(0,P,Q), U, V2)
        #y[,,j] <- rmatrixt(n=1,df,mean=matrix(0,P,Q), U2, V2)
      }
      
    }
    
    
    ###### ##### ##### ##### ##### ##### 
    ###### ##### ##### ##### ##### #####
    
    x.vec <- matrix(0,N1,P*Q)
    y.vec <- matrix(0,N2,P*Q)
    for (j in 1:N1){  x.vec[j,]<- matrix(x[,,j], 1, P*Q) }
    for (j in 1:N2){  y.vec[j,]<- matrix(y[,,j], 1, P*Q) }
    xy.vec<- rbind(x.vec, y.vec)
    d <- dist(xy.vec)
    E<-mstree(d, 3)
    sigma1 = sqrt(0.5*median(d^2))
    
    
    #gamma =0.1  #0.01
    gamma =1/sigma1
    pvalue.Detc<- twosample.test(x,y,gamma,B)
    
    
    ###### ##### MMD##### ##### #####
    lab  <- c(rep(1,N1), rep(2,N2))  
    K_gau= dnorm(as.matrix(dist(xy.vec,diag=T,upper=T)), mean=0, sd=sigma1)
    pvalue.MMD<- mmd2test(K_gau, lab)$p.value
    
    ###### ##### Energy##### ##### #####
    pvalue.Energy= eqdist.etest(d, sizes=c(N1, N2), distance=TRUE, R = B)$p.value
    
    ###### ##### Graph based Test##### ##### #####
    pvalue.Gtest<-(g.tests(E,1:N1,(N1+1):(N1+N2))$weighted)$pval.approx
    
    ###### ##### Robust Graph based##### ##### #####
    pvalue.rgtest=rg.test(x.vec,y.vec, n1=N1, n2=N2, k = 5, 
                          weigh.fun=weiMax, perm.num = 0)$asy.wei.pval
    ###### ##### Ball##### ##### #####
    pvalue.ball= bd.test(x = x.vec, y = y.vec)$p.value
    
    ###### ##### Chen and Qin##### ##### #####
    pvalue.CQ=apval_Chen2010(x.vec,y.vec)$pval

    
    pvalue.Kim=Per_Test_Kim2020Aos(x.vec,y.vec,99)
    ###### ##### ##### ##### ##### ##### 
 
    temp<-data.frame(Detc=pvalue.Detc,  MMD=pvalue.MMD, Energy=pvalue.Energy,
                     Graph=pvalue.Gtest, RGtest=pvalue.rgtest, Ball=pvalue.ball,
                     CQ=pvalue.CQ,Kim=pvalue.Kim)
    pvalues<-rbind(pvalues, temp)
    cat("Inner loop",loop,"Out loop",i,'\n')
  }
  
  
  
  apha=0.05
  p_Detc=sum(pvalues$Detc<apha)/nsim
  p_MMD=sum(pvalues$MMD<apha)/nsim
  p_Energy=sum(pvalues$Energy<apha)/nsim
  p_Graph=sum(pvalues$Graph<apha)/nsim
  p_RGtest=sum(pvalues$RGtest<apha)/nsim
  p_Ball=sum(pvalues$Ball<apha)/nsim
  p_CQ=sum(pvalues$CQ<apha)/nsim
  p_Kim=sum(pvalues$Kim<apha)/nsim
  
  pp=data.frame(Detc=p_Detc,MMD=p_MMD, Energy=p_Detc,
                Graph=p_Graph,RGtest=p_RGtest, Ball=p_Ball,
                CQ=p_CQ,Kim=p_Kim)
  power<-rbind(power, pp)
  N_seq[i]=rho2
  
}


res=cbind(N_seq,power)
print(res)

par(mfrow = c(1, 1))  
matplot(N_seq, power,type="b",cex.lab=1.2, lty = 2:6, xlab=expression(u[2]),
        ylab="empirical power",lwd=3.5, cex = 1.5, pch = 21:25, col = c(1,2,4,6,7) )


abline(h = apha,  col = "lightgray") 
axis(4,  at=apha,  col.axis = "red", lwd =2)
legend("topleft", c("DPCD", "Energy","MMD_gaus","MMD_lap","Ball"), lty = 2:6,lwd=3.5,
       cex = 1.0, pch = 21:25,col = c(1,2,4,6,7) )