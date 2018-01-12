
diabetes <- read.csv("diabetes.csv")

## BLASSO

library(monomvn) ## install.packages("monomvn")

diabetes_design <- model.matrix(prog ~ . - 1, data = diabetes)
diabetes_design_std <- scale(diabetes_design, center = TRUE, scale = TRUE)

diabetes_response <- diabetes$prog - mean(diabetes$prog)

diabetes_blasso <- blasso(X = diabetes_design_std, y = diabetes_response, T = 500, thin = 200)

boxplot(diabetes_blasso$beta, horizontal = TRUE, las = 1)
abline(v=0, lty = 1)

hist(rowSums(abs(diabetes_blasso$beta)), probability = TRUE)

## BOOTSTRAP

diabetes_n <- nrow(diabetes)
bootstrap_idx <- replicate(1000, sample(diabetes_n, size = diabetes_n, replace = TRUE), simplify = FALSE)

library(glmnet)
diabetes_cv <- cv.glmnet(x = diabetes_design_std, y = diabetes_response)
lambda_cv <- cv.glmnet(x = diabetes_design_std, y = diabetes_response)$lambda.min

diabetes_boot <- sapply(bootstrap_idx, function(idx) 
  as.matrix(coef(glmnet(x = diabetes_design_std[idx,], y = diabetes_response[idx], lambda = lambda_cv))))

boxplot(t(diabetes_boot[-1,]), horizontal = TRUE, las = 1)
abline(v=0, lty = 1)


