#EX 1
#' # Picloram data
#' 
#' Picloram is a herbicide used for plant control and Turner
#' et. al. (1992) examined the effect of varying doses of Picloram for
#' control of tall larkspur (ridderspore). Four different doses of
#' Picloram (0, 1.1, 2.2, 4.5 kg/ha) were used on in total 313 plants
#' located in three areas (called `replicate` in the dataset). For
#' each of the plants it was recorded whether the plant was killed or
#' not.

picloram <- read.table("data/picloram.txt", header=TRUE)
picloram
require(ggplot2)
qplot(dose, dead/total, data=picloram, colour=replicate, geom=c("point","line"))

#' The overall questions are: How lethal is the herbicide at different
#' doses and is there a difference in effect of dosis in the different
#' locations.
#' 
#' What would a good model be? Fit such a model? Discuss pros and cons
#' of different models. Does it fit well to data? Interpret the
#' parameters. What do you conclude?

#' NB: A sigmoid curve can look like this
x <- seq(0.01, 0.99, 0.01)
qplot(log(x/(1-x)), x, xlab="Logit(theta)", ylab=expression(theta), geom="line")

#' These code fragments may come in handy.
picloram <- transform(picloram, elogit=log((dead + 0.5)/(total - dead + 0.5)))
qplot(dose, elogit, data=picloram, colour=replicate, geom=c("point","line"))

logreg1 <- glm( cbind( dead, total - dead ) ~ replicate + dose, family=binomial, 
                data=picloram)

coef(summary(logreg1))



