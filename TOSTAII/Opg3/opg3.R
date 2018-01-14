data(Orthodont, package="nlme")
ort <- Orthodont
head(ort)
levels(ort$Subject)
library(ggplot2)
library(tidyverse)
qplot(age, distance, group=Subject, data=ort, color=Subject) + geom_path() + facet_grid(~Sex)

osum <- ort %>% group_by(Sex, age) %>% summarize(m=mean(distance), v=var(distance))
qplot(m, v, data=osum)

library(geepack)
geeglm(formula = distance~age,data = ort,id=Subject,corstr="unstructured")
unique(ort$age)
lm(formula = distance~age,data = ort)
unique(ort$Subject)

ort1 <- ort[order(ort$Subject,ort$age),]

fit_normal_gee <- function(formula, id, corstr="independence", phi=NULL,data, w=NULL){
  #data_temp <- data[order(id)]
  if(!is.null(phi)){
    if((length(phi)>1|phi<0)[1]){
      stop("phi should be just one possitive number")
    }
  }
  deparsed_id <- deparse(substitute(id))
  y <- model.response(model.frame(formula, data=data))
  X <- model.matrix(formula, data=data)
  subjects <- unique(data[deparsed_id])
  n_subs <- length(subjects)
  n_vars <- ncol(X)
  n_obs <- length(y)
  n_times <- n_obs/n_subs
  print(n_times)
  #print(subjects)
  if (is.null(w))
    w <- rep(1, length(y))
  conv.eps <- 1e-12
  ## Get started; intitial value of z and beta
  lm. <- lm.fit(x = X, y = y)
  beta_star <- lm.$coefficients
  print(beta_star)
  k <- 0
  repeat{
    pearsons_resid <- (y-X%*%beta_star)/1#Variansen er sat til 1
    #print(pearsons_resid)
    phi <- t(pearsons_resid)%*%pearsons_resid/(n_obs-n_vars)
    print(phi)
    #print(pearsons_resid)
    rmatrix <- diag(1,4)
    #1/((n_subs-n_vars)*phi)
    #for(i in subjects){
    for(i in 0:(n_subs-1)){
      individ_res <- pearsons_resid[(1:n_times)+(i*n_times)]
      #individ_res <- pearsons_resid[subjects==i]
      #print(X[subjects==i,])
      corr <- individ_res%*%t(individ_res)
      #diag(corr) <- 0
      rmatrix <- rmatrix+corr
      #print(individ_res)
      #print(corr-t(corr))
    }
    rmatrix <- (1/((n_subs-n_vars)*phi[1,1]))*rmatrix
    diag(rmatrix) <- 1
    #print(rmatrix)
    #v_i <- phi*rmatrix
    print(rmatrix)
    #return(rmatrix)
    rbig <- kronecker(diag(1,n_subs),rmatrix)
    mu_star <- X%*%beta_star
    beta_hat <- beta_star+solve(t(X)%*%solve(rbig)%*%X)%*%t(X)%*%solve(rbig)%*%(y-mu_star)
    if(sum(abs(beta_star-beta_hat)<conv.eps)==ncol(X)){break}
    beta_star <- beta_hat
    print(beta_hat)
    k <- k+1
    print(k)
    if(k>100){break}
  }
  
}
fit_normal_gee(formula = distance~age*Sex,data = ort,id=Subject,phi = 1)
kk <- fit_normal_gee(formula = distance~age,data = ort,id=Subject,phi = 1)


for (j in 1:nobsses) {
  for (k in 1:nobsses) {
    if(k!=j){
      corr[i,k]<- individ_res[j]*individ_res[k]
    }
  }
}
    repeat{
      pearson_resid <- (y-mu_hat)/sqrt(V(mu_hat)/w)
      if(sum(abs(beta_star-beta_hat)<conv.eps)==ncol(X)){break}
    }
}
  
  ## Iterate until beta no longer changes
  repeat{
    eta_star <- X%*%beta_star
    mu_star <- g.inv(eta_star)
    z_star <- g(mu_star)+g.der(mu_star)*(y-mu_star)
    v_star <- (g.der(mu_star)^2*V(mu_star))/w
    #if(link=="identity"){H <- diag(rep(1,length(y)))}else{H <- diag(g.der(mu_star)[,1])}
    sigma <- diag(v_star[,1])
    #sigma <- diag(v_star)
    #r <- H%*%(y-mu_star)
    #z=X%*%beta_star+r
    beta_hat <- solve(t(X)%*%solve(sigma)%*%X)%*%t(X)%*%solve(sigma)%*%z_star
    if(sum(abs(beta_star-beta_hat)<conv.eps)==ncol(X)){break}
    beta_star <- beta_hat
    # break the loop when beta stops changing between succesive
    # iterations.
  }
  ## Compute mu and working weights v after final iteration
  mu_hat <- g.inv(X%*%beta_hat)
  ## Estimate phi if necessary
  p <- ncol(X)
  if(!is.null(phi)){
    phi_hat <- phi
  }else{
    phi_hat <- (1/(nrow(X)-p))*sum(((y-mu_hat)^2)/(V(mu_hat)/w))
  }
  ## Compute estimated covariance matrix of regression parameters
  vcov <- phi_hat*solve(t(X)%*%solve(sigma)%*%X)
  ## Compute Pearson residuals
  pearson_resid <- (y-mu_hat)/sqrt(V(mu_hat)/w)
  
  
  out <- list( 
    coef = beta_hat,
    vcov = vcov,
    phi  = phi_hat,
    resid= pearson_resid,
    fit  = mu_hat,
    p    = p)
  out 
}
multiplot(qplot(log2(conc), 1/lot1, data=ct),
          qplot(log2(conc), log(lot1), data=ct),
          qplot(log2(conc), lot1, data=ct),
          cols=3)

## OPG 3
g1a <- glm(lot1 ~ log2(conc), family=Gamma("inverse"), data=ct)
fga <- fit_gamma(lot1 ~ log2(conc),link = "inverse",data=ct)
g1a$coefficients-fga$coef
summa1 <- summary(g1a)
summa1$cov.scaled-fga$vcov
resid(g1a,type = "pearson")
fga$resid

g1b <- glm(lot1 ~ log2(conc), family=Gamma("log"), data=ct)
fgb <- fit_gamma(lot1 ~ log2(conc),link = "log",data=ct)
g1b$coefficients
fgb$coef

g1c <- glm(lot1 ~ log2(conc), family=Gamma("identity"), data=ct)
fgc <- fit_gamma(lot1 ~ log2(conc),link = "identity",data=ct)
g1c$coefficients
fgc$coef


## OPG 4
plot(fga$fit,fga$resid)
plot(log2(ct$conc), 1/ct$lot1)
lines(1/fga$fit,col="red")
plot(fgb$fit,fgb$resid)
plot(log2(ct$conc), log(ct$lot1))
lines(log(fgb$fit),col="red")
plot(fgc$fit,fgc$resid)
plot(log2(ct$conc), ct$lot1)
lines(fgb$fit,col="red")

##OPG 5
library(Matrix)
diag(diag(1,4),matrix(0,nrow = 4,ncol = 4))
a <- diag(1,4)
b <- matrix(0,nrow = 4,ncol = 4)

b[2,2] <- diag(1,2)


g1d <- glm(lot1 ~ log2(conc)+I(log2(conc)^2), family=Gamma("inverse"), data=ct)
fgd <- fit_gamma(lot1 ~ log2(conc)+I(log2(conc)^2),link = "inverse",data=ct)

summa2 <- summary(g1d)

wtest <- (g1d$coef^2)/diag(summa2$cov.scaled)
pchisq(wtest,lower.tail = FALSE,df = 3)
wtest <- (fgd$coef^2)/diag(fgd$vcov)
pchisq(wtest,lower.tail = FALSE,df = 3)

wtest <- (fga$coef^2)/diag(fga$vcov)
pchisq(wtest,lower.tail = FALSE,df = 2)