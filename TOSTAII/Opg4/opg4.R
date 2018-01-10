
therm <- read.csv("thermal.csv", header=T)


fit_rational_nls <- function(y=1,x=2,p=2,q=2, data){
  data_temp <- cbind(data[1:(1+p+q),y],data[1:(1+p+q),x])
  colnames(data_temp) <- c("y","x")
  data_temp <- as.data.frame(data_temp)
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
  print(nls_temp)
}
fit_rational_nls(data=therm)
fit_rational_nls(y="tec",x="temp",data=therm,p=4,q=5)


