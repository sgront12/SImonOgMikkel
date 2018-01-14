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

### OPG 2
fit_normal_gee <- function(formula, id, corstr="independence", phi=NULL,data, w=NULL){
  #data_temp <- data[order(id)]
  
  if(!is.null(phi)){
    set_phi <- TRUE
    if((length(phi)>1|phi<0)[1]){
      stop("phi should be just one possitive number")
    }
  }else{set_phi <-  FALSE}
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
  if(corstr=="unstructured"){
      repeat{
      pearsons_resid <- (y-X%*%beta)/1#Variansen er sat til 1
      if(set_phi == FALSE){
        phi <- drop(t(pearsons_resid)%*%pearsons_resid/(n_obs-n_vars))
      }
      

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
    cov_y <- (y-X%*%beta_hat)%*%t(y-X%*%beta_hat)
    for(i in 1:(n_subs-1)){
      cov_y[((n_subobs*i)+1):n_obs,i:(n_subobs*i)] <- 0
      cov_y[i:(n_subobs*i),((n_subobs*i)+1):n_obs] <- 0
    }
    I_1 <- 1/phi^2*t(X)%*%rbig_inv%*%cov_y%*%rbig_inv%*%X
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
  I_0 <- solve(t(X)%*%X)
  cov_y <- (y-X%*%beta)%*%t(y-X%*%beta)
  for(i in 1:(n_subs-1)){
    cov_y[((n_subobs*i)+1):n_obs,i:(n_subobs*i)] <- 0
    cov_y[i:(n_subobs*i),((n_subobs*i)+1):n_obs] <- 0
  }
  I_1 <- 1/phi^2*t(X)%*%cov_y%*%X
  return(list("coef" = beta,
              "phi" = phi,
              "p" = length(beta),
              "rmatrix"=diag(1,n_subobs),
              "Sigma_m" = I_0,
              "vcov" = I_0%*%I_1%*%I_0, #Sigma_e,
              "fit" = X%*%beta,
              "resid" = y-X%*%beta))
  
}
### OPG 3
fit1 <- fit_normal_gee(formula = distance~age+Sex,data = ort,id=Subject,phi = 1,corstr="Independent")
fit3 <- fit_normal_gee(formula = distance~age+Sex,data = ort,id=Subject,corstr="unstructured")
fit1$phi
fit3$phi
fit1$coef
fit3$coef
sqrt(diag(fit3$vcov))


### OPG 4
#Parallele
par(mfrow=c(1,2))
fit4 <- fit_normal_gee(formula = distance~age+Sex,data = ort,id=Subject,corstr="unstructured")
plot(ort$age,ort$distance,col=as.numeric(ort$Sex)+4)
## male line
abline(fit4$coef["(Intercept)",],fit4$coef["age",],col="blue")
## female line
abline(fit4$coef["(Intercept)",]+fit4$coef["SexFemale",],fit4$coef["age",],col="red")
fit5 <- fit_normal_gee(formula = distance~age*Sex,data = ort,id=Subject,corstr="unstructured")
plot(ort$age,ort$distance,col=as.numeric(ort$Sex)+4)
## male line
abline(fit5$coef["(Intercept)",],fit5$coef["age",],col="blue")
## female line
abline(fit5$coef["(Intercept)",]+fit5$coef["SexFemale",],fit5$coef["age",]+fit5$coef["age:SexFemale",],col="red")

### OPG 5
plot(fit4$fit,fit4$resid)
plot(fit5$fit,fit5$resid)

### OPG 6
wtest <- (fit5$coef^2)/diag(fit5$vcov)
pchisq(wtest,lower.tail = FALSE,df = 4)

### OPG 7
fit2 <- geeglm(formula = distance~age*Sex,data = ort,id=Subject,corstr="unstructured")
fit2$geese$vbeta
fit2$geese$
fit2
fit1$coef
fit1$phi





