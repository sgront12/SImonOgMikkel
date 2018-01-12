library(c060)

set.seed (3)
x <- matrix(rnorm (80 * 500 ,0 ,1) ,80 ,500)
y <- x[1:80 ,1:500] %*% c(rep (3 ,2),rep(-3,3),rep (.1 ,495))

res <- stabpath(y,x)
sel <- stabsel(res ,error =0.05 , type="pfer",pi_thr =0.6)
sel

plot(res)

## Diabetes

diabetes_stabpath <- stabpath(y = diabetes_response, x = diabetes_design_std)
diabetes_stable_selected <- stabsel(diabetes_stabpath ,error =0.05 , type="pfer",pi_thr =0.6)
diabetes_stable_selected

plot(diabetes_stabpath)
