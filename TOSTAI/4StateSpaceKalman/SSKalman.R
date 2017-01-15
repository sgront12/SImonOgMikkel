# Kalman-filter

D = c(4,5,4,1,0,4,3,4,0,6,
      3,3,4,0,2,6,3,3,5,4,5,3,1,4,4,1,5,5,3,4,2,5,2,2,3,4,
      2,1,3,2,1,1,1,1,1,3,0,0,1,0,1,1,0,0,3,1,0,3,2,2,0,1,
      1,1,0,1,0,1,0,0,0,2,1,0,0,0,1,1,0,2,2,3,1,1,2,1,1,1,
      1,2,4,2,0,0,0,1,4,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0)
yr = 1851:1962
N = length(D)
points(yr, D,type="p")

#1. Judging from looking at data, is classical ARIMA modelling a viable road to model these data? Why? Why not?


#2. Is there any evidence of serial correlation in data?


#3. Using the dlm package define a local level model (also known as random walk plus noise model).

var(D)
mean(D)
ms <- dlm(m0=4, C0=10, FF=1, V=var(D), GG=1, W=10)
flt <- dlmFilter(y, mod=ms)
plot(yr, y)
lines(yr, dropFirst(flt$m))

#4. Plot the filtered values; plot the forecasts.


#5. Fit the variance parameters to data.


#6. Plot the filtered values; plot the forecasts again.


#7. Looking at data, there is some evidence of a change point. Fit a local level model that allows for such a change point.


#8. Plot the filtered values; plot the forecasts again.


#9. What do you conclude from all this?