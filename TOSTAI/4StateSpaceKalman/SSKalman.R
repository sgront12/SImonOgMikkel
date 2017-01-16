# Kalman-filter

library(dlm)
library(magrittr)
library(ggplot2)
library(dplyr)
library(broom)
library(changepoint)


D = c(4,5,4,1,0,4,3,4,0,6,
      3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,
      2,1,3,2,1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,
      1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,
      1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0)
yr = 1851:1962
N = length(D)
plot(yr, D,type="p")

#1. Judging from looking at data, is classical ARIMA modelling a viable road to model these data? Why? Why not?
#Nej arima kan ikke håndtere changepoints

#2. Is there any evidence of serial correlation in data?
#Det ser sådan ud, men det føles ikke sådan

#3. Using the dlm package define a local level model (also known as random walk plus noise model).

var(D)
ms <- dlm(m0=4, C0=10, FF=1, V=2000, GG=1, W=10)
flt <- dlmFilter(D, mod=ms)


#4. Plot the filtered values; plot the forecasts.
yrp10<-1851:1972
Dp10<-append(D,c(rep(NA,10)))
plot(yrp10, Dp10)
lines(yr, dropFirst(flt$m))
?dlmForecast
for1<-dlmForecast(flt, nAhead = 10, method = c("plain", "svd"), sampleNew = FALSE)
for1fp10<-append(rep(NA,length(D)),for1$f)
lines(yrp10,for1fp10,col=c("red"))


#5. Fit the variance parameters to data.
build <- function(parm){
  spec <- dlm(m0=4, C0=10, FF=1, V=exp(parm[1]), GG=1, W=exp(parm[2]))
  spec
}
build(c(1, 1))

p <- dlmMLE(D, parm=c(1, 1), build)$par
ms.fit <- build( p )
build2 <- function(parm){
  spec <- dlm(m0=4, C0=exp(parm[3]), FF=1, V=exp(parm[1]), GG=1, W=exp(parm[2]))
  spec
}
build2(c(1, 1, 1))

p2 <- dlmMLE(D, parm=c(1, 1, 1), build2)$par
ms2.fit <- build2( p2 )
#6. Plot the filtered values; plot the forecasts again.
par(mfrow=c(3,1))
plot(yr, D)
dlmFilter(D, ms)$m %>% dropFirst %>% lines(yr, .)
plot(yr, D)
dlmFilter(D, ms.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, D)
dlmFilter(D, ms2.fit)$m %>% dropFirst %>% lines(yr, .)

#7. Looking at data, there is some evidence of a change point. Fit a local level model that allows for such a change point.
cpt.var(D)
cpt.mean(D)
cpt.meanvar(D)
x <- rep(1, length(D))
x[30:45] <- 5
x
plot(yr, x)

X <- matrix(x)
ms  <- dlm(m0=4, C0=10, FF=1, V=2000, GG=1, JW=1, W=10, X=X)

flt <- dlmFilter(D, mod=ms)
plot(yr, y)
lines(yr, dropFirst(flt$m))
abline(v=1886, col="red", lwd=2)

build3 <- function(parm){
  Z <- X
  Z[,1] <- exp(parm[1])*X[,1]
  spec <- dlm(m0=4, C0=10, FF=1, V=exp(parm[2]), GG=1, JW=1, W=10, X=Z)
  spec
}
build3(1)

p3 <- dlmMLE(D, parm=c(1,1), build3)$par
ms3.fit <- build3( p3 )
plot(yr, D)
dlmFilter(D, ms3.fit)$m %>% dropFirst %>% lines(yr, .)




#8. Plot the filtered values; plot the forecasts again.
par(mfrow=c(2,2))
plot(yr, D)
dlmFilter(D, ms)$m %>% dropFirst %>% lines(yr, .)
plot(yr, D)
dlmFilter(D, ms.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, D)
dlmFilter(D, ms2.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, D)
dlmFilter(D, ms3.fit)$m %>% dropFirst %>% lines(yr, .)

#9. What do you conclude from all this?