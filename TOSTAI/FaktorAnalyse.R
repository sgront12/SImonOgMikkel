#Faktor analyse
library(corrplot)
#Funktion, som estimere parametrene i en FA model vha. EM-algoritmen
faEM=function(y,nfac,tol=1e-4,stand=F){
  #Data centreres(og standardiseres, hvis stand==T)
  x=as.matrix(scale(y,scale=stand))
 
  #Begyndelsesparametrene bestemmes
  C=cov(x)
  eig_val=eigen(C)$values[1:nfac]
  eig_vec=eigen(C)$vectors[,1:nfac]
  Lambda_t1=eig_vec%*%diag(sqrt(eig_val),nfac)
  Psi_t1=diag(diag(C-Lambda_t1%*%t(Lambda_t1)))
  
  #Matricer til Lambda_t og Psi_t
  Lambda_t=matrix(0,ncol(x),nfac)
  Psi_t=matrix(0,ncol(x),ncol(x))
  
  #Vektorer til E-trinet
  mean_f=matrix(0,nfac,nrow(x))
  mean_ff=matrix(0,nfac,nfac)
  
  #Til M-trinet
  xMeanf=matrix(0,ncol(x),nfac)
  xx=t(x)%*%x
  
  #Vektor til log-likelihood
  ll=c()
  
  #Variabel til bestemmelse af antal iterationer
  k=0
  
  #EM-algoritmen
  while(all(abs(Lambda_t1-Lambda_t) > tol) || all(abs(Psi_t1-Psi_t) > tol)){
    mean_ff[,]=0
    Lambda_t=Lambda_t1
    Psi_t=Psi_t1
    var_f=solve(diag(1,nfac)+t(Lambda_t)%*%solve(Psi_t, Lambda_t))
    inv_var_x=solve(Lambda_t %*% t(Lambda_t)+Psi_t)
    
    #Log-likelihood
    ll[k+1] = (nrow(x)/2)*log(det(inv_var_x))-(1/2)*sum(diag(xx %*%inv_var_x))
    
    
    #E-trin
    for(i in 1:nrow(x)){
      mean_f[,i]=t(Lambda_t)%*%inv_var_x%*%x[i,]
      mean_ff=mean_ff+var_f+mean_f[,i]%*%t(mean_f[,i])
    }
    
    #M-trin
    xMeanf=xx%*%inv_var_x%*%Lambda_t
    
    Lambda_t1=xMeanf%*%solve(mean_ff)
    Psi_t1=diag(diag(xx-Lambda_t1%*%t(xMeanf))/nrow(x))
    
    k=k+1
  }
  
  #Den sidste log-likelihood
  inv_var_x=solve(Lambda_t1 %*% t(Lambda_t1)+Psi_t1)
  ll[k+1]=(nrow(x)/2)*log(det(inv_var_x))-(1/2)*sum(diag(t(x) %*% x %*%inv_var_x))
  
  #De estimerede parametre returneres
  return(list(Lambda=Lambda_t1,Psi=Psi_t1,LogLike=ll,iterationer=k))
}

#FA model med 1 faktor
X=read.csv("~/Desktop/math.csv") #Datasaet
stand=T #Data standardiseres
tol=1.e-6 #Tolerance
nfac=1 #Antal faktorer

FF1=faEM(X,nfac,tol,stand=stand)

FF1$Lambda; FF1$Psi; FF1$iterationer

#FA model med 2 faktorer
nfac=2 #Antal faktorer

FF2=faEM(X,nfac,tol,stand=stand)

FF2$Lambda; FF2$Psi; FF2$iterationer

varimax(FF2$Lambda)$loadings[,1:2]
corrplot(cor(X), order = "hclust", tl.col='black', tl.cex=1.5)
