#' ## Dynamic linear models examples
#' ### based on the `dlm` package

library(dlm)
library(magrittr)
library(ggplot2)
library(dplyr)
library(broom)

#' ## Nile data
#' 
#' Annual flow of river Nile at Aswan in 1871-1970, in 10^8 m^3, “with apparent
#' changepoint near 1898”
#' 
#' ### Random walk + noise aka Local level model
#'
#' ### Nile data

y <- as.numeric(Nile)
yr  <- 1871:1970
plot(yr, y)

#' Create dlm model object; parameters are pretty much arbitrary
ms <- dlm(m0=1100, C0=10, FF=1, V=1000, GG=1, W=10)
flt <- dlmFilter(y, mod=ms)
plot(yr, y)
lines(yr, dropFirst(flt$m))

#' Get ready to estimate parameters from data
build <- function(parm){
  spec <- dlm(m0=1100, C0=10, FF=1, V=exp(parm[1]), GG=1, W=exp(parm[2]))
  spec
}
build(c(1, 1))

p <- dlmMLE(y, parm=c(1, 1), build)$par
ms.fit <- build( p )

#' More elaborate: Estimate 3 parameters
build2 <- function(parm){
  spec <- dlm(m0=1100, C0=exp(parm[3]), FF=1, V=exp(parm[1]), GG=1, W=exp(parm[2]))
  spec
}
build2(c(1, 1, 1))

p2 <- dlmMLE(y, parm=c(1, 1, 1), build2)$par
ms2.fit <- build2( p2 )

par(mfrow=c(1,3))
plot(yr, y)
dlmFilter(y, ms)$m %>% dropFirst %>% lines(yr, .)
plot(yr, y)
dlmFilter(y, ms.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, y)
dlmFilter(y, ms2.fit)$m %>% dropFirst %>% lines(yr, .)


#' Nile data - with change point around year 1900
#' 
#' Increase variance by a factor >1 around year 1900 to allow for level change:
#' 
x <- rep(1, length(y))
x[20:40] <- 5
x
plot(yr, x)

X <- matrix(x)
ms  <- dlm(m0=1100, C0=10, FF=1, V=1000, GG=1, JW=1, W=10, X=X)

flt <- dlmFilter(y, mod=ms)
plot(yr, y)
lines(yr, dropFirst(flt$m))
abline(v=1898, col="red", lwd=2)

build3 <- function(parm){
  Z <- X
  Z[,1] <- exp(parm[1])*X[,1]
  spec <- dlm(m0=1100, C0=10, FF=1, V=exp(parm[2]), GG=1, JW=1, W=10, X=Z)
  spec
}
build3(1)

p3 <- dlmMLE(y, parm=c(1,1), build3)$par
ms3.fit <- build3( p3 )
plot(yr, y)
dlmFilter(y, ms3.fit)$m %>% dropFirst %>% lines(yr, .)

par(mfrow=c(2,2))
plot(yr, y)
dlmFilter(y, ms)$m %>% dropFirst %>% lines(yr, .)
plot(yr, y)
dlmFilter(y, ms.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, y)
dlmFilter(y, ms2.fit)$m %>% dropFirst %>% lines(yr, .)
plot(yr, y)
dlmFilter(y, ms3.fit)$m %>% dropFirst %>% lines(yr, .)



#' ## Pig growth
#'
#' ### Linear (local) growth model
data(dietox, package="doBy")
head(dietox)

#' Look at subset of data; approximately linear growth 
diet <- subset(dietox, Evit=="Evit000" & Cu=="Cu000",
               select=c("Pig", "Weight", "Time"))
ggplot(diet, aes(x=Time, y=Weight, color=Pig)) + geom_point() + geom_line()

#' But not quite linear growth;
#' Grow rate increases with about .2 kg per week:
diet %>%
  group_by(Pig)  %>%
  mutate(dWeight = Weight - lag(Weight)) -> diet.diff
ggplot(diet.diff,
       aes(x=Time, y=dWeight, color=Pig)) + geom_point() + geom_line()
mm <- lm(dWeight ~ Time, data=diet.diff)
mm %>% tidy
mm %>% glance

#' Bottom line: Approximately linear growth but regression coefficients evolve
#' over time (one possible model)

#' Pick one pig to look at
y <- subset(diet, Pig=="4601")$Weight
plot(y)

#' Specify a model; no date at this point
ms <- dlm(m0=c(20,5), C0=diag(1,2),
          FF=matrix(c(1,0), nr=1), V=0.01,
          GG=matrix(c(1, 0, 1, 1), nr=2), W=diag(c(.1, .1)))

#' Combine model and data in filtering
flt <- dlmFilter(y, mod=ms)
par(mfrow=c(1,2))
plot(y)
lines(dropFirst(flt$m[,1]))
plot(dropFirst(flt$m[,2]))

#' Get ready to fit (some) parameters 
build <- function(parm){
  spec <- dlm(m0=c(20, 5), C0=diag(1,2),
              FF=matrix(c(1, 0), nr=1), V=exp(parm[1]),
              GG=matrix(c(1, 0, 1, 1), nr=2),
              W=diag(c(exp(parm[2]), exp(parm[3]))))
  spec
}
#' Create a model object
build( c(0, 0, 0) )

#' Fit model to data and create model object
p <- dlmMLE(y, parm=c(0,0,0), build)$par
ms.fit <- build(p)

#' Compare model fits
flt2 <- dlmFilter(y, mod=ms.fit)
par(mfrow=c(2,2))
plot(y)
lines(dropFirst(flt$m[,1]))
plot(dropFirst(flt$m[,2]))
plot(y)
lines(dropFirst(flt2$m[,1]))
plot(dropFirst(flt2$m[,2]))

#' saveGIF(
#' for (i in 1:12){
#'     tt <- 1:15
#'     use <- 1:i
#'     f <- dlmFilter(y[use], mod=ms.fit)
#'     filt <- dropFirst(f$m[,1])
#'     fore <- dlmForecast(f, nAhead=15)$f  %>% as.numeric
#'     fore <- fore[1:(length(tt)-length(use))]
#'     comb <- c(filt, fore )
#'     pch <- 1+c(rep(1, length(filt)), rep(2, length(fore)))
#'     plot(comb, pch=pch, ylim=c(10, 140))
#'     points(y, col="red", lwd=2)
#'     lines(filt)
#' },
#' movie.name="pig-growth2.gif")
#' 















ms <- dlm(m0=c(0,0), C0=diag(1,2),
          FF=matrix(c(1,1), nr=1), V=1.4,
          GG=diag(1,2), W=diag(0,2))

f <- dlmFilter(y, mod=ms)
plot(y)
lines(dropFirst(f$m[,1]))

build <- function(parm){
  ms <- dlm(m0=c(1200,-2.7), C0=diag(1000,2),
            FF=matrix(c(parm[1], parm[2]), nr=1), V=exp(parm[3]),
            GG=diag(1,2), W=diag(0,2))
  ms
}
mm <- build(c(1,1,0))

dlmMLE(y, parm=c(1, 2, 0), mm)$par -> p
p
m <- build( p )
f <- dlmFilter(y, m)
pred <- f$m %*% p 
head(pred)
plot(y)
lines(dropFirst(pred))



dlm(FF=matrix(c(1,0), nr=1),
    V=1.4,
    GG=matrix(c(1, 0, 1, 1), nr=2),
    W = diag(c(0, 0.2)),
    m0=rep(0, 2),
    C0=10*diag(2))


N <- 100
y <- cumsum(rnorm(N))
par(mfrow=c(1,2))
pacf(y)
acf(y)
arima(y)

y=c(4,5,4,1,0,4,3,4,0,6,
    3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,
    2,1,3,2,1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,
    1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,
    1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0)
yr=1851:1962


ms <- dlm(m0=0, C0=10, FF=1, V=1.4, GG=1, W=.1)
is.dlm(ms)

f <- dlmFilter(y, mod=ms)
plot(y)
lines(dropFirst(f$m))

dlmMLE(ms)

build <- function(parm){
  ms <- dlm(m0=0, C0=10, FF=1, V=exp(parm[1]), GG=1, W=exp(parm[2]))
  ms
}

m <- build(c(1.2, 0.2))

dlmMLE(y, parm=c(1, 2), m)$par -> p
m <- build( p )

plot(y)
dlmFilter(y, m)$m  %>% dropFirst  %>% lines
lines(dropFirst(f$m))



#' Go to wide format:
dietw <-
  reshape(diet, direction="wide",
          idvar="Time",
          timevar="Pig",
          v.names = "Weight")

dw <- apply(dietw, 2, diff) %>% as.data.frame
dw$Time <- dietw$Time[-1]

matplot(dw[,-1], type="l")

dw$Time <- factor(dw$Time)
dl <- melt(dw)
dl$Time <- as.numeric(dl$Time)

dw <- apply(dietw, 2, diff) %>% as.data.frame
dw$Time <- dietw$Time[-1]

dl <- reshape(dw, direction="long",
              varying=2:ncol(dw),
              timevar="Pig")
dl <- transform(dl, Pig=factor(Pig))
head(dl, 13)              

