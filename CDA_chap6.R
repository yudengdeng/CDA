#chap6a
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
dta$revcolor <-  factor(dta$color, levels = c("5", "2", "3", "4")) 
dta$revspine <-  factor(dta$spine, levels = c("3", "1", "2")) 
dta$weight <- dta$weight/1000

library(MASS)
n = dim(dta)[1]
#step(model_full,direction = "backward") #AIC
model.final <- step(model_full, k=log(n),direction = "backward",trace = TRUE) #BIC
