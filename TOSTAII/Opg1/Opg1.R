ct <- read.table("data/clottingTime.txt", header=T)
ct
library(ggplot2)
source("multiplot.R")
fit_gamma <- function(formula, data, link="inverse", phi=NULL, w=NULL){
  y <- model.response(model.frame(formula, data=data))
  X <- model.matrix(formula, data=data)
  f <- Gamma(link)
  g <- f$linkfun
  g.inv <- f$linkinv
  g.der <- Deriv::Deriv(g, "mu")
  V <- f$variance
  if (is.null(w))
    w <- rep(1, length(y))
  conv.eps <- 1e-12
  ## Get started; intitial value of z and beta
  lm. <- lm.fit(x = X, y = g(y))
  
  beta_star <- lm.$coefficients
  
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
  phi_hat <- (1/(nrow(X)-p))*sum(((y-mu_hat)^2)/(V(mu_hat)/w))
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
multiplot(qplot(log2(conc), lot1, data=ct),
          qplot(log2(conc), 1/lot1, data=ct),
          qplot(log2(conc), log(lot1), data=ct), cols=3)

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

fga$resid
fgb$resid
fgc$resid




a <- c(1,2,3,4,5)
b <- c(2,3,4,77,6)
a-b
if(sum(abs(a-b))<2){print("tis")}


if(require("microbenchmark")) {
  mb <- microbenchmark(lm(y~X), lm.fit(X,y), .lm.fit(X,y))
  print(mb)
