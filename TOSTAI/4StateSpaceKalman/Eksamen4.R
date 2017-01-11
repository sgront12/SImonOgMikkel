D=c(4,5,4,1,0,4,3,4,0,6,3,3,4,0,2,6,3,3,5,4,5,3,
    1,4,4,1,5,5,3,4,2,5,2,2,3,4,2,1,3,2,1,1,1,1,
    1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,1,1,0,1,
    0,1,0,0,0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,
    1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0)
yr = 1851:1962
N  =length(D)
plot(yr, D)

#1. Judging from looking at data, is classical ARIMA modelling a viable road
#to model these data? Why? Why not?

#Svar: Virker umiddelbart poissonfordelt, hvilket taler imod ARIMA, men eftersom
###### vi ikke kender antallet af kulminer, kunne dette tale for ARIMA.

#2. Is there any evidence of serial correlation in data?
plot(yr[2:1112]-yr[1])

library(dlm)
?dlm
ms <- dlm(m0=2, C0=1, FF=1, V=1, GG=2, W=1)
fil = dlmFilter(D, mod=ms)
plot(yr, D)
lines(yr, dropFirst(fil$m))
length(D)
