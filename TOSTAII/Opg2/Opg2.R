# Ex 2 Penalised regression
leukemia <- read.csv("leukemia.csv", header=T)
library(c060)
library(glmnet)

x <- as.matrix(leukemia[,-(1:2)])
fit1 <- glmnet(x,leukemia$class,family = "binomial",alpha = 1)
summary(fit1)
plot(fit1)
fit1$dim
head(x)
cv.lasso <- cv.glmnet(x, leukemia$class, family='binomial', alpha=1, parallel=FALSE, standardize=TRUE)
head(coef(cv.lasso))
plot(cv.lasso)
plot(cv.lasso$glmnet.fit, xvar="lambda", label=TRUE)
cv.lasso$lambda.min
cv.lasso$lambda.1se
fittet1 <- cv.lasso$glmnet.fit
fittet1$lambda[2]
leukemiaX <- leukemia[,-(1:2)]
colnames(leukemiaX[which(coef(cv.lasso, s=fittet1$lambda[2])[-1,]!=0)])





## BOOTSTRAP
leukemia_design <- model.matrix(class ~ . - 1-X, data = leukemia)
leukemia_design_std <- scale(leukemia_design, center = TRUE, scale = TRUE)

leukemia_response <- leukemia$class
leukemia_n <- nrow(leukemia)
bootstrap_idx <- replicate(1000, sample(leukemia_n, size = leukemia_n, replace = TRUE), simplify = FALSE)

leukemia_cv <- cv.lasso
  #cv.glmnet(x = leukemia_design_std, y = leukemia_response)
lambda_cv <- cv.lasso$lambda.min
  #cv.glmnet(x = leukemia_design_std, y = leukemia_response)$lambda.min

head(coef(cv.lasso))

leukemia_boot <- sapply(bootstrap_idx, function(idx) 
  as.matrix(coef(glmnet(x = leukemia_design_std[idx,], y = leukemia_response[idx], family='binomial', standardize=TRUE, lambda = lambda_cv))))

#boxplot(t(leukemia_boot[-1,]), horizontal = TRUE, las = 1)
#abline(v=0, lty = 1)
rows <- rowSums(leukemia_boot[-1,]!=0)
hist(rows)
rows <- rows[rows!=0]
hist(rows)
rows <- rows[rows>=200]
hist(rows)
sort(unique(rowSums(leukemia_boot[-1,]!=0)))
 sum(rowSums(leukemia_boot[-1,]!=0)>=700)
leukemiaX <- leukemia[,-(1:2)]
colnames(leukemiaX[,(rowSums(leukemia_boot[-1,]!=0)>=700)])
##Stability selection

leukemia_stabpath <- stabpath(y = leukemia_response, x = leukemia_design,family="binomial")
leukemia_stable_selected <- stabsel(leukemia_stabpath ,error =0.05 , type="pfer",pi_thr =0.6)
leukemia_stable_selected$stable

plot(leukemia_stabpath)

#cv.glm
colnames(leukemiaX[which(coef(cv.lasso, s=fittet1$lambda[2])[-1,]!=0)])
#bootstrap
colnames(leukemiaX[,(rowSums(leukemia_boot[-1,]!=0)>=700)])
#stability selection
leukemia_stable_selected$stable

