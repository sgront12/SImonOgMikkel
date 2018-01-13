data(Orthodont, package="nlme")
ort <- Orthodont
library(geepack)
ort1 <- ort[order(rnorm(108)),]
ort2 <- ort[1:16,]
ort3 <- ort2[order(rnorm(16)),]
ort4 <- ort3[-2,]

geeglm(formula = distance~age,data = ort,id=Subject,corstr="unstructured")

fit_normal_gee <- function(formula, id,mtime, corstr="independence", phi=NULL,data, w=NULL){
  #data_temp <- data[order(id)]
  if(!is.null(phi)){
    if((length(phi)>1|phi<0)[1]){
      stop("phi should be just one possitive number")
    }
  }
  y <- model.response(model.frame(formula, data=data))
  n_obs <- length(y)
  X <- model.matrix(formula, data=data)
  n_vars <- ncol(X)
  deparsed_id <- deparse(substitute(id))
  subjects <- unique(data[deparsed_id])
  n_subs <- length(subjects)
  deparsed_mtime <- deparse(substitute(mtime))
  times <- unique(data[deparsed_mtime])
  n_times <- length(times)
  otimes <- times[order(times)]
  n_zerorows <- 0
  #print(times)
  if(corstr=="unstructured"){
    X_ordered <- matrix(0,nrow = 0,ncol = n_vars)
    y_ordered <- matrix(0,nrow = 0,ncol = 1)
    for(i in subjects){
      X_i <- matrix(0,nrow = n_times,ncol = n_vars)
      y_i <- matrix(0,nrow = n_times,ncol = 1)
      for (j in 1:n_times) {
        who <- (data[deparsed_id]==i&data[deparsed_mtime]==otimes[j])
        #print(who)
        #print(length(X[who,]))
        if(length(X[who,])!=0){
          #print(X[who,])
        X_i[j,] <- X[who,]
        y_i[j,] <- y[who]
        }else(n_zerorows <- n_zerorows+1)
      }
      X_ordered <- rbind(X_ordered,X_i)
      y_ordered <- rbind(y_ordered,y_i)
      #print(X_ordered)
      #print(y_ordered)
      
    }
  }
  #print(subjects)
  if (is.null(w))
    w <- rep(1, length(y))
  conv.eps <- 1e-12
  ## Get started; intitial value of z and beta
  lm. <- lm.fit(x = X_ordered, y = y_ordered)
  beta_star <- as.matrix(drop(lm.$coefficients))
  rownames(beta_star) <- colnames(X)
  print(beta_star)
  k <- 0
  repeat{
    pearsons_resid <- (y_ordered-X_ordered%*%beta_star)/1#Variansen er sat til 1
    #print(pearsons_resid)
    phi <- (t(pearsons_resid)%*%pearsons_resid)/(n_obs-n_vars)
    print(phi)
    #print(pearsons_resid)
    rmatrix <- diag(0,n_times)
    #1/((n_subs-n_vars)*phi)
    for(i in 0:(n_subs-1)){
      individ_res <- pearsons_resid[(1:n_times)+(i*n_times)]
      #print(X_ordered[(1:n_times)+(i*n_times),])
      #print(individ_res)
      corr <- individ_res%*%t(individ_res)
      #corr <- cov2cor(corr)
      #diag(corr) <- 0
      rmatrix <- rmatrix+corr
      #print(rmatrix)
      #print(individ_res)
      #print(corr)
    }
    #print(rmatrix)
    rmatrix <- (1/((n_subs-n_vars)*phi[1,1]))*rmatrix
    #print(rmatrix)
    #rmatrix <- cov2cor(rmatrix)
    diag(rmatrix) <- 1
    
    #v_i <- phi*rmatrix
    #print(rmatrix)
    #return(rmatrix)
    rbig <- kronecker(diag(1,n_subs),rmatrix)
    #print(rbig)
    mu_star <- X_ordered%*%beta_star
    beta_hat <- beta_star+solve(t(X_ordered)%*%solve(rbig)%*%X_ordered)%*%t(X_ordered)%*%solve(rbig)%*%(y_ordered-mu_star)
    if(sum(abs(beta_star-beta_hat)<conv.eps)==ncol(X)){break}
    beta_star <- beta_hat
    #print(beta_hat)
    k <- k+1
    print(k)
    if(k>200){break}
  }
  print(rmatrix)
  return(beta_hat)
}
fit_normal_gee(formula = distance~age,data = ort,id=Subject,mtime = age,corstr="unstructured")
kk <- fit_normal_gee(formula = distance~age,data = ort,id=Subject,phi = 1)

