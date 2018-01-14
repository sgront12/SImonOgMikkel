library(ggplot2)
library(tidyverse)
library(geepack)

data(Orthodont, package="nlme")
ort <- Orthodont
head(ort)
levels(ort$Subject)

qplot(age, distance, group=Subject, data=ort, color=Subject) + geom_path() + facet_grid(~Sex)
osum <- ort %>% group_by(Sex, age) %>% summarize(m=mean(distance), v=var(distance))
qplot(m, v, data=osum)


fit_normal_gee <- function(formula, id, corstr="independence", phi=NULL,data, w=NULL){
  #data_temp <- data[order(id)]
  if(!is.null(phi)){
    if((length(phi)>1|phi<0)[1]){
      stop("phi should be just one possitive number")
    }
  }
  deparsed_id <- deparse(substitute(id))
  y <- model.response(model.frame(formula, data=ort))
  X <- model.matrix(formula, data=data)
  subjects <- unique(data[deparsed_id])
  n_subs <- length(subjects)
  n_vars <- ncol(X)
  n_obs <- length(y)
  n_subobs <- n_obs/n_subs
  conv.eps <- 1e-14
  ## Get started; intitial value of z and beta
  lm. <- lm.fit(x = X, y = y)
  beta <- lm.$coefficients
  
  #Beviset for, at det er en lortefunktion, der blev brugt tidligere.
  #v = c()
  #for(i in 1:n_subs){
  #  v <- c(v,as.vector(match(pearsons_resid[subjects==subjects[i]], pearsons_resid)))
  #}
  
  repeat{
    pearsons_resid <- (y-X%*%beta)/1#Variansen er sat til 1
    
    phi <- drop(t(pearsons_resid)%*%pearsons_resid/(n_obs-n_vars))

    rmatrix <- diag(1,n_subobs)
    
    for(i in 1:n_subs){
      individ_res <- pearsons_resid[(n_subobs*(i-1)+1):(i*n_subobs)]
      corr <- individ_res%*%t(individ_res)
      rmatrix <- rmatrix+corr
    }
    rmatrix <- 1/((n_subs-n_vars)*phi)*rmatrix
    diag(rmatrix) <- 1
    rbig_inv <- solve(kronecker(diag(1,n_subs),rmatrix))
    
    beta_hat <- beta+solve(t(X)%*%rbig_inv%*%X)%*%t(X)%*%rbig_inv%*%(y-X%*%beta)
    if(sum(abs(beta-beta_hat)<conv.eps)==ncol(X)){break}
    beta <- beta_hat
  }
  #Needs some fixing, those values are off the charts
  I_0 <- solve(t(X)%*%rbig_inv%*%X)
  I_1 <- 1/phi^2*t(X)%*%rbig_inv%*%(y-X%*%beta_hat)%*%t(y-X%*%beta_hat)%*%rbig_inv%*%X
  Sigma_e <- I_0%*%I_1%*%I_0

  return(list("coef" = beta_hat,
              "phi" = phi,
              "p" = length(beta_hat),
              "rmatrix"=rmatrix,
              "Sigma_m" = I_0,
              "vcov" = Sigma_e,
              "fit" = X%*%beta_hat,
              "resid" = y-X%*%beta_hat))
}
fit1 <- fit_normal_gee(formula = distance~age+Sex,data = ort,id=Subject,phi = 1)


fit2 <- geeglm(formula = distance~age+Sex,data = ort,id=Subject,corstr="unstructured")
fit2$geese$vbeta
