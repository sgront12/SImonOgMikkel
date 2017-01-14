# source("http://bioconductor.org/biocLite.R");
# biocLite(c("graph","RBGL","Rgraphviz"))
# install.packages("gRbase", dependencies=TRUE)
# install.packages("gRain", dependencies=TRUE)
# install.packages("gRim", dependencies=TRUE)

library(gRbase)
library(gRain)
library(gRim)
library(Rgraphviz)


##Opgave 1
dg1<-dag(~S+L|S+X|L:S+B|S+D|L:B)
dg2<-dag(~S|X+L|S+X|L+B|S+D|L:B)
par(mfrow=c(1,2))
plot(dg1)
plot(dg2)
b1<-as(dg1, "matrix")
b2<-as(dg2, "matrix")

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

algoritme(b1)
algoritme(b2)

##Opgave 2
dg<-dag(~smoke+lung|smoke+xray|lung:smoke+bronc|smoke+dysp|bronc:lung)
plot(dg)
data(chestSim1000,package="gRbase")
head(chestSim1000)

#1. Extract the necessary CPTs from data, and construct the Bayesian network.
extractCPT(chestSim1000,dg,smooth = 0)
#Frekvenstabel over data
ftab=xtabs(~smoke+lung+bronc+dysp+xray,data=chestSim1000)
#Test af betinget uafhaengighed
# P(s,x,l,b,d)
ciTest(ftab,set=c("dysp","xray","smoke","lung","bronc"))
# P(D|L,B,S)
ciTest(ftab,set=c("dysp","smoke","lung","bronc")) # vi kan tage smoke ud
ciTest(ftab,set=c("dysp","lung","bronc")) # vi kan ikke tage lung ud
ciTest(ftab,set=c("dysp","bronc","lung")) # vi kan ikke tage bronc ud
# P(X|L,D,S)
ciTest(ftab,set=c("xray","smoke","lung","bronc"))# vi kan tage smoke ud
ciTest(ftab,set=c("xray","lung","bronc")) # vi kan ikke tage lung ud
ciTest(ftab,set=c("xray","bronc","lung")) # vi kan tage bronc ud
# P(L,D,S)
ciTest(ftab,set=c("lung","bronc","smoke")) # vi kan tage bronc ud
# P(L|S)
ciTest(ftab,set=c("lung","smoke")) # vi kan ikke tage smoke ud
# P(L|S)
ciTest(ftab,set=c("bronc","smoke")) # vi kan ikke tage smoke ud
# vi ender med # P(D|L,B,S)
cpt<-compileCPT(extractCPT(chestSim1000,dg,smooth = 0))
bn<-grain(cpt)
plot(bn)
#2. What does information about dysp tell us about smoke, i.e. what is the conditional distri- bution smoke given dysp?
querygrain(bn, nodes=c("smoke","dysp"), type="conditional")
#3. If we know smoke, what does additional information about bronc tell us about lung? That is, what is the conditional distribution of lung given smoke, and what is the conditional distribution of lung given smoke and bronc?
querygrain(bn, nodes=c("lung","smoke"), type="conditional")
querygrain(bn, nodes=c("lung","smoke","bronc"), type="conditional")
#4. If we know smoke and dysp, what does additional information about bronc tell us about lung?

#5. Sketch the message passing algorithm (CollectEvidence and DistributeEvidence) for this specific example.




