# Selvom han ikke anbefaler load, kunne jeg ikke få den ind med attatch
load("carcass_fa.RData")


# Jeg starter med et korellationsplot
library(corrplot)
corrplot(cor(carcass_fa), order = "hclust", tl.col='black', tl.cex=.75)
# Og hvad er det så vi kan læse ud fra det her, at de klumpede firkanter symboliserer, lad os sige 3 grupper

# Så standardiseres der (her skal uddybes hvorfor)
carcass <- as.data.frame(scale(carcass_fa))

# Vi pruger psych pakken, der står godt nok den er beregnet til pprykiatriske tests, men lad os se hvad vi får ud af det
load(psych)
nfact_carcass <- psych::nfactors(carcass,rotate = "none")
# Jeg forstår ikke alle de her plots, en det ser ud til i den "Emperical BIC" at 3 faktorer er sagen

# det her er vist en anden måde at bestemme antallet af faktorer på
plot(nscree_carcass <- nFactors::nScree(carcass))
# igen ser det ud ti at "optimal coordinates" er 3
nscree_carcass

# Så finder vi faktorerne
FA_none <- factanal(carcass, factors = 3, rotation = "none", na.action = na.omit)
FA_none$loadings

# Her er jeg ikke helt med
sum(FA_none$loadings[,1]^2)/32 ## S is correlation matrix, hence: tr(S) = 1*p = 32


# Vi ser på nogle egenværdier
(carcass_eigen <- eigen(cov(carcass))$value)
# Og akkurart som egenværdierne var større end 1 i 6 tilfælde med self_rate, er de det i 3 tilfælde med carcass

# Her er jeg ikke lige med
FA_none$uniquenesses
loadings_distant = FA_none$loadings[1,]
(communality_distant = sum(loadings_distant^2))
(uniqueness_distant = 1-communality_distant)


# Men nu vil vi måske rotere
# Vi starter med lige at plotte
par(mfrow=c(1,1))
plot(FA_none$loadings[,1:2], type="n") # set up plot
text(FA_none$loadings[,1:2],labels=rownames(FA_none$loadings),cex=.7)


# Så roterer vi hen til "varimax"
FA_varimax = factanal(carcass, factors = 3, rotation = "varimax",na.action = na.omit, scores = "regression")

plot(FA_varimax$loadings[,1:2], type="n") # set up plot
text(FA_varimax$loadings[,1:2],labels=rownames(FA_varimax$loadings),cex=.7)


# Der er blevet en mere tydelig gruppering, men "weight" ser mere ensom ud

pairs(FA_varimax$loadings[,1:3],panel = function(x,y) text(x,y, labels=rownames(FA_varimax$loadings)))



# til sidst den der oblimax
library(GPArotation)
(FA_oblimax = factanal(self_rate, factors = 3, rotation = "oblimax"))
# Det er godt nok en lav p-værdi, men det er den også med kun 2 faktorer

