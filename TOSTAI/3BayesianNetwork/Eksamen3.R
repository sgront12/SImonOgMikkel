source("http://bioconductor.org/biocLite.R");
biocLite(c("graph","RBGL","Rgraphviz"))
install.packages("gRbase", dependencies=TRUE)
install.packages("gRain", dependencies=TRUE)
install.packages("gRim", dependencies=TRUE)

library(gRbase)
library(gRain)
library(gRim)
library(Rgraphviz)

dg1            <-            dag(~S +  L|S+ X|L:S +  B|S  +    D|L:B)
dg2            <-            dag(~S|X + L|S + X|L +  B|S  +    D|L:B)
par(mfrow=c(1,2))
plot(dg1)
plot(dg2)
b<-as(dg2, "matrix")
b
nrow(b)
rowSums(b)

algoritme <- function(x){ #x er en matrix
  while(length(x)>1){
    sum = rowSums(x)
    if(any(sum==0)==TRUE){
      for(i in 1:nrow(x)){
        if(sum[i]==0){ 
          x <- x[-i,-i]
        }
      }
    }else{
      return(print("Det er ikke DAG!"))
    }
  }
  if(x[1]==0 || length(x)==0){
    print("Det er en DAG!")
  }else{print("Det er ikke en DAG!")}
}

algoritme(b)
b
