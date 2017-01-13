mnEM=function(x, tol){
  #x er en nx2 matrix bestaaende af realiseringer af en multinomiel fordeling.
  #tol er tolerancen. 
  k=0;
  #Antal realiseringer
  n=dim(x)[1]
  #Gør vektoren (n_00,n_10,n_01,n_11,n_0+,n_1+,n_+0,n_+1) klar.
  #Data clean-up
  counts = c(numeric(8))
  i=1
  while(i <= n){
    if(is.na(x[i,1])){x[i,1]=2}
    if(is.na(x[i,2])){x[i,2]=2}
    if(x[i,1]==0){
      if(x[i,2]==0){ counts[1]<- counts[1]+1;}
      else if(x[i,2]==1){ counts[3]<- counts[3]+1;}
      else{ counts[5]<- counts[5]+1;}
    }
    else if(x[i,1]==1){
      if(x[i,2]==0){ counts[2]<- counts[2]+1;}
      else if(x[i,2]==1){ counts[4]<- counts[4]+1;}
      else{ counts[6]<- counts[6]+1;}
    }
    else{
      if(x[i,2]==0){ counts[7]<- counts[7]+1;}
      else if(x[i,2]==1){ counts[8]<- counts[8]+1;}
    }
    i <- i+1;
  }
  
  #E og M vektorer
  #Vektor til pi_ij
  pi_t1=c(1/4,1/4,1/4,1/4)
  pi_t=numeric(4)
  
  #Vektor til den sufficiente stikproevefunktion
  # s_t = (nI_00, nI_01, nI_10, nI_11, NJ_00, NJ_01, NJ_10, NJ_11)
  s_t=numeric(8)
  
  while(all((abs(pi_t1-pi_t)) > tol) || k>1000){
    pi_t = pi_t1
    #E-trin
    s_t[1] = counts[5]*pi_t[1]/(pi_t[1]+pi_t[2])
    s_t[2] = counts[5]*pi_t[2]/(pi_t[1]+pi_t[2])
    s_t[3] = counts[6]*pi_t[3]/(pi_t[3]+pi_t[4])
    s_t[4] = counts[6]*pi_t[4]/(pi_t[3]+pi_t[4])
    
    s_t[5] = counts[7]*pi_t[1]/(pi_t[1]+pi_t[3])
    s_t[6] = counts[8]*pi_t[2]/(pi_t[2]+pi_t[4])
    s_t[7] = counts[7]*pi_t[3]/(pi_t[1]+pi_t[3])
    s_t[8] = counts[8]*pi_t[4]/(pi_t[2]+pi_t[4])
    
    #M-trin
    pi_t1[1] = (counts[1]+s_t[1]+s_t[5])/sum(counts)
    pi_t1[2] = (counts[2]+s_t[2]+s_t[6])/sum(counts)
    pi_t1[3] = (counts[3]+s_t[3]+s_t[7])/sum(counts)
    pi_t1[4] = (counts[4]+s_t[4]+s_t[8])/sum(counts)
    
    k <- k+1;
  }
  return(list(tal = counts, Estep = s_t, itt = k, pi = pi_t1));
}

#Anvendelse af ovenstaaende funktion
library(readr)
IJdata <- read_csv("C:/Users/Mikke/Desktop/SImonOgMikkel/trunk/TOSTAI/1EMalgorithm/IJdata.csv")
View(IJdata)
tol=1.e-10 #Tolerance
X=IJdata #Datasaet

EM1=mnEM(X,tol)
EM1$Estep
sum(EM1$pi)
EM1$pi
EM1$itt