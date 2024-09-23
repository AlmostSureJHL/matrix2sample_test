#===============================================
# Generate Covariance Matrix
#===============================================
generate_sig = function(n,rho){
  if(rho==0)
    sig <- diag(n)     
  else{
    a <- seq(1:n)
    b <- matrix(a,n,n)
    c <- t(b)
    q <- abs(b-c)
    tmp <- q*log(rho)
    sig <- exp(tmp)
  }       
  return(sig)
}

generate_band = function(n,rho){
  sig<-diag(1,n,n)
  if(rho==0)
    sig <- diag(n)     
  else{
    for(i in 1:n){
      for(j in 1:n){
        if(abs(i-j)==1){
          sig[i,j]<-rho
        }
        if(abs(i-j)==2){
          sig[i,j]=rho/2
        }
      }
    }
  }       
  return(sig)
}


#===============================================
# two sample test  
#===============================================
twosample.test<- function(x,y,gamma,B){
  n<- dim(x)[3]
  m<- dim(y)[3]
  p<- dim(x)[1]
  q<- dim(x)[2]
  z<- rep(0,p*q*(n+m))
  dim(z) <- c(p,q,n+m)
  z[,,1:n]<-x
  z[,,(n+1):(n+m)]<-y
  K_mat<- matrix(NA,n+m,n+m)
  One<- as.matrix(c(rep(1,n)/n,-rep(1,m)/m))
  L<- One%*%t(One)

    I_q<- diag(q)
    for (i in 1:(n+m)){
      for (j in 1:(n+m)){
        zz_ij=matrix(z[,,i]-z[,,j], p, q)
        K_mat[i,j]<-(det(gamma*t(zz_ij)%*%zz_ij+I_q))^(-1/2)
      }
    }
    
  

  P_mat<- tr(K_mat%*%L)
  K_sampM<- matrix(NA,n+m,n+m)
  
  samp <- lapply(1:B, function(x) sample(1:(n+m)))
  reps <- sapply(1:B,
                 function(v) {
                   K_sampM <- K_mat[samp[[v]],samp[[v]]]
                   res<- tr(K_sampM%*%L)
                  return(res)
                 })
  pval <- (1 + length(which(reps>P_mat))) / (1 + B)
}



twosample.sym.test<- function(x,y,gamma,B){
  n<- dim(x)[3]
  m<- dim(y)[3]
  p<- dim(x)[1]
  q<- dim(x)[2]
  z<- rep(0,p*q*(n+m))
  dim(z) <- c(p,q,n+m)
  z[,,1:n]<-x
  z[,,(n+1):(n+m)]<-y
  K_mat<- matrix(NA,n+m,n+m)
  One<- as.matrix(c(rep(1,n)/n,-rep(1,m)/m))
  L<- One%*%t(One)

    I_q<- diag(q)
    for (i in 1:(n+m)){
      for (j in 1:(n+m)){
        zz_ij=matrix(z[,,i]-z[,,j], p, q)
 
       DET= prod(eigen(I_q-2*gamma*zz_ij*1i, only.values=TRUE)$values)
        K_mat[i,j]<-DET^(-1/2)
      }
    }
    
 



  P_mat<- Re( tr(K_mat%*%L))
  ###print(  tr(K_mat%*%L) )
  K_sampM<- matrix(NA,n+m,n+m)
  
  samp <- lapply(1:B, function(x) sample(1:(n+m)))
  reps <- sapply(1:B,
                 function(v) {
                   K_sampM <- K_mat[samp[[v]],samp[[v]]]
                   res<-Re(  tr(K_sampM%*%L) )
                  return(res)
                 })
  pval <- (1 + length(which(reps>P_mat))) / (1 + B)
}


MSingular<- function(x,y){
  n<- dim(x)[3]
  m<- dim(y)[3]
  p<- dim(x)[1]
  q<- dim(x)[2]
  z<- rep(0,p*q*(n+m))
  dim(z) <- c(p,q,n+m)
  z[,,1:n]<-x
  z[,,(n+1):(n+m)]<-y
  MaxSingular<- matrix(NA,n+m,n+m)

  for (i in 1:(n+m)){
    for (j in 1:(n+m)){
      zz_ij<- matrix(z[,,i]-z[,,j], p, q)
      MaxSingular[i,j]<- max(svd(zz_ij)$d )
    }
  }
  MaxSingular_1<- matrix(MaxSingular,1,(n+m)^2)
  M<- median(MaxSingular_1)
  return(M)
}
 



Kim2020Aos<-function(x.vec,y.vec){
  n<-nrow(x.vec)
  m<-nrow(y.vec)
  Second<-0
  Third<-0
  for(i in 1:m){
    Temp1y<-x.vec-y.vec[i,]
    Temp2y<-Temp1y%*%t(Temp1y)
    N_Temp1y<-sqrt(diag(Temp2y))
    N_Temp1y_M<-N_Temp1y%*%t(N_Temp1y)
    Temp3y<-Temp2y/N_Temp1y_M
    diag(Temp3y)<-1
    Second<-Second+sum(acos(Temp3y)/(n^2*m))
    }
  for(i in 1:n){
    Temp1x<-y.vec-x.vec[i,]
    Temp2x<-Temp1x%*%t(Temp1x)
    N_Temp1x<-sqrt(diag(Temp2x))
    N_Temp1x_M<-N_Temp1x%*%t(N_Temp1x)
    Temp3x<-Temp2x/N_Temp1x_M
    diag(Temp3x)<-1
    Third<-Third+sum(acos(Temp3x)/(m^2*n)) 
  }
    
  Stat<-1/3-0.5*pi*Second-0.5*pi*Third

  return(Stat)
}

Per_Test_Kim2020Aos<-function(x.vec,y.vec,B){
  n<-nrow(x.vec)
  m<-nrow(y.vec)
  xy.vec<-rbind(x.vec,y.vec)
  samp <- lapply(1:B, function(x) sample(1:(n+m)))
  Stat<-Kim2020Aos(x.vec,y.vec)
  reps <- sapply(1:B,
                 function(v) {
                   Ord<-samp[[v]]
                   xSamp_ind<-Ord[1:n]
                   ySamp_ind<-Ord[-(1:n)]
                   
                   x.vec_sample<-xy.vec[xSamp_ind,]
                   y.vec_sample<-xy.vec[ySamp_ind,]
                   
                   res<-Kim2020Aos(x.vec_sample,y.vec_sample)
                   return(res)
                 })
  pval <- (1 + length(which(reps>Stat))) / (1 + B)
  
}

