ct <- read.table("data/clottingTime.txt", header=T)
ct
library(ggplot2)
source("multiplot.R")
multiplot(qplot(log2(conc), lot1, data=ct),
          qplot(log2(conc), 1/lot1, data=ct),
          qplot(log2(conc), log(lot1), data=ct), cols=3)

g1a <- glm(lot1 ~ log2(conc), family=Gamma("inverse"), data=ct)
g1a

g1b <- glm(lot1 ~ log2(conc), family=Gamma("log"), data=ct)
g1b

g1c <- glm(lot1 ~ log2(conc), family=Gamma("identity"), data=ct)
g1c

f<- Gamma("inverse")
summary(f)
?lm.fit
require(utils)

formula <- lot1 ~ log2(conc)
y <- model.response(model.frame(formula, data=ct))
X <- model.matrix(formula, data=ct)
y
lm. <- lm.fit (x = X, y = y)

set.seed(129)



if(require("microbenchmark")) {
  mb <- microbenchmark(lm(y~X), lm.fit(X,y), .lm.fit(X,y))
  print(mb)



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
  lm. <- lm.fit(x = X, y = y)
g(y)
  
  beta_star <- lm.$coefficients

  ## Iterate until beta no longer changes
  repeat{
    eta_star <- X%*%beta_star
    mu_star <- g.inv(eta_star)
    z_star <- g(mu_star)+g.der(mu_star)*(y-mu_star)
    v_star <- (g.der(mu_star)^2*V(mu_star))/w
    if(link=="identity"){H <- diag(rep(1,length(y)))}
    else{H <- diag(g.der(mu_star)[,1])}
    sigma <- diag(v_star[,1])
    #r <- H%*%(y-mu_star)
    #z=X%*%beta_star+r
    beta_hat <- solve(t(X)%*%solve(sigma)%*%X)%*%t(X)%*%solve(sigma)%*%z_star
    if(sum(abs(beta_star-beta_hat)<conv.eps)==ncol(X)){break}
    beta_star <- beta_hat
    # break the loop when beta stops changing between succesive
    # iterations.
  }
    ## Compute mu and working weights v after final iteration
    ## Estimate phi if necessary
    ## Compute estimated covariance matrix of regression parameters
    ## Compute Pearson residuals
    out <- list( ## Fill in your values
      coef = NULL,
      vcov = NULL,
      phi  = NULL,
      resid= NULL,
      fit  = NULL,
      p    = NULL)
  out 
}

a <- c(1,2,3,4,5)
b <- c(2,3,4,77,6)
a-b
if(sum(abs(a-b))<2){print("tis")}
