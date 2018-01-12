
therm <- read.csv("thermal.csv", header=T)


fit_rational_nls <- function(y=1,x=2,p=2,q=2, data){
  data_temp <- cbind(data[1:(1+p+q),y],data[1:(1+p+q),x])
  colnames(data_temp) <- c("y","x")
  matrix_temp <- as.matrix(rep(1,(1+p+q)))
  names_temp <- c("alpha0")
  y_name <- colnames(data[y])
  #print(y_name)
  x_name <- colnames(data[x])
  #print(x_name)
  formula_temp <- paste(y_name,"~alpha0",sep = "")
  for (i in 1:p) {
    matrix_temp <- cbind(matrix_temp,data_temp[,2]^i)
    names_temp <- append(names_temp,paste("alpha",i,sep = ""))
    formula_temp <- paste(formula_temp,"+alpha",i,"*",x_name,"^",i,sep="")
  }
  for (i in 1:q) {
    matrix_temp <- cbind(matrix_temp,-data_temp[,1]*data_temp[,2]^i)
    names_temp <- append(names_temp,paste("beta",i,sep = ""))
    formula_temp <- paste(formula_temp,"-beta",i,"*",y_name,"*",x_name,"^",i,sep="")
  }
  #print(matrix_temp)
  startval <- solve(matrix_temp,data_temp[,1])
  names(startval) <-names_temp
  #print(startval)
  #print(formula_temp)
  formula_temp <- as.formula(formula_temp)
  #print(formula_temp)
   nls_temp <- nls(formula_temp
                  , data=data
                  , start = startval)
  return(nls_temp)
}
fit1 <- fit_rational_nls(data=therm)
summary(fit1)
fit2 <- fit_rational_nls(y="tec",x="temp",data=therm,p=2,q=2)
summary(fit2)

kk <- 1
AIC <- rep(NA,25)
BIC <- rep(NA,25)
name_round <- rep(NA,25)
for (k in 1:100) {
  for (j in 1:100) {
    name_round[kk] <- paste(k,j)
    #print(name_round[kk])
    try(fit2 <- fit_rational_nls(y="tec",x="temp",data=therm,p=k,q=j),silent = TRUE)
    #if((k==3 &j==5)|(k==4 &j==4)|(k==5 &j==5)){print("noGo")}
    #else{fit2 <- fit_rational_nls(y="tec",x="temp",data=therm,p=k,q=j)
    AIC[kk] <- AIC(fit2)
    BIC[kk] <- BIC(fit2)
    #print(BIC(fit2))
    #}
    kk <- kk+1
  }
  
}
paste("AIC siger",name_round[min(which(min(AIC, na.rm = TRUE)==AIC))])
paste("BIC siger",name_round[min(which(min(BIC, na.rm = TRUE)==BIC))])

fit3 <- fit_rational_nls(y="tec",x="temp",data=therm,p=2,q=6)
fit4 <- fit_rational_nls(y="tec",x="temp",data=therm,p=3,q=6)
plot(therm,pch = 19,col="black")
points(therm$temp,residuals(fit2)+therm$tec,col="green")
points(therm$temp,residuals(fit3)+therm$tec,col="blue")
points(therm$temp,residuals(fit4)+therm$tec,col="red")