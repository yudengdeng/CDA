#chap 4a

snoring   <- c(0,2,4,5)
disease <- c(24,35,21,30)
total <- c(1355,603,192,224)
result <- glm(cbind(disease,total)~ snoring, family = binomial(link ="identity"))
summary(result)


result <- glm(cbind(disease,total)~ snoring, family = binomial(link ="logit"))
summary(result)

#logisticCurv <- function(x){exp(a+b*x) / (1+exp(a+b*x))}
#a = -2
#b =0.4
#curve(logisticCurv, col = 1)
#curve(function2, col = 2, add = TRUE)
#curve(function3, col = 3, add = TRUE)
#curve(function4, col = 4, add = TRUE)

#chap 4b
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
w <- dta$width
wsq <- dta$width^2
sa <- dta$satell

result <- glm(sa~ w, family = poisson(link ="log"))
summary(result)

result <- glm(sa~ w, family = poisson(link ="identity"),start=c(0.5,0.5), maxit=500)
summary(result)

result <- glm(sa~ w + wsq, family = poisson(link ="log"))
summary(result)

#install.packages("MASS")
library(MASS)
result <- glm.nb(sa ~ w, link="log")
summary(result)
result$theta #shape parameter, i.e. k

result <- glm.nb(sa ~ w, link="identity",start=c(0.5,0.5), maxit=500)
result$theta #shape parameter, i.e. k
summary(result)
