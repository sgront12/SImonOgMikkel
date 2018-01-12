library(relaxo) ## install.packages("relaxo")

## Assumes bayesian_bootstrap run before this

diabetes_relaxo <- relaxo(Y=diabetes_response, X = diabetes_design_std)

plot(diabetes_relaxo)

diabetes_relaxo_cv <- cvrelaxo(Y=diabetes_response, X = diabetes_design_std)
print(diabetes_relaxo_cv$phi)
print(diabetes_relaxo_cv$lambda)

fitted_values_relaxo <- predict(diabetes_relaxo_cv)
fitted_values_lasso <- predict(diabetes_cv, newx = diabetes_design_std)[,1]
plot(fitted_values_relaxo,diabetes_response, ylab = "Observed", xlab = "Predicted")
points(fitted_values_lasso,diabetes_response, col = 3)
legend("topleft", bty = "n", col = c(1,3), c("Relaxed LASSO", "LASSO"), pch = 1)
abline(c(0,1))

var(fitted_values_relaxo-diabetes_response)
var(fitted_values_lasso-diabetes_response)
