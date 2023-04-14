
#chap 5a
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
y <- ifelse(dta$satell > 0, 1, 0) # y = a binary indicator of satellites
w <- dta$width
result <- glm(y ~ w, family=binomial(link=logit))
summary(result)

#group the data

dta$y <- ifelse(dta$satell > 0, 1, 0) # y = a binary indicator of satellites
dta$w = NA
dta[dta$width > 20.25,]$w = 20.75
dta[dta$width > 21.25,]$w = 21.75
dta[dta$width > 22.25,]$w = 22.75
dta[dta$width > 23.25,]$w = 23.75
dta[dta$width > 24.25,]$w = 24.75
dta[dta$width > 25.25,]$w = 25.75
dta[dta$width > 26.25,]$w = 26.75
dta[dta$width > 27.25,]$w = 27.75
dta[dta$width > 28.25,]$w = 27.75
dta[dta$width > 29.25,]$w = 28.75

dta$cnt <- 1
dta_agg <- aggregate(dta,by = list(dta$w),FUN = sum)
dta_agg$logit = log((dta_agg$y+0.5)/(dta_agg$cnt-dta_agg$y+0.5))
names(dta_agg)[1] = "wNew"
plot(dta_agg$wNew,dta_agg$logit)


#plot 95% CIs
lower = 20
upper = 33
confint(result)
dtaW_disp <- seq(lower,upper,0.1) 
newdata <- data.frame(w = dtaW_disp)
crab.predict <- predict(result,newdata, type="response", se=T)

plot(dtaW_disp,crab.predict$fit, axes=F, type="l", xlim=c(lower,upper),
     ylab="Probability of satellite", xlab="") 
axis(2, at=seq(0,1,.2))
axis(1, at=seq(20,32,2))
lines(dtaW_disp,crab.predict$fit-1.96*crab.predict$se.fit,lty=3)
lines(dtaW_disp,crab.predict$fit+1.96*crab.predict$se.fit,lty=3)



#Goodness-of-fit for Logistic Regression for Ungrouped Data
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
dta$w <- dta$width
dta$y <- ifelse(dta$satell > 0, 1, 0) # y = a binary indicator of satellites
#First, we create a table with the successes and failures per width category 
cont.table<-xtabs(~w+y, data=dta)
w.unique <-as.numeric(attr(cont.table,"dimnames")$w)
matrix.succ.fail<-structure(.Data=cont.table,dim=c(66,2))[,2:1] #observed successes and failures

w.cut <- cut(w.unique, breaks=c(0,seq(23.25, 29.25),Inf), left.include=T)
observed <- apply(matrix.succ.fail,2,aggregate,by=list(W=w.cut),sum)
observed <- matrix(c(observed[[1]][,ncol(observed[[1]])],
                     observed[[2]][,ncol(observed[[2]])]), ncol = 2)

fit.logit <- glm(matrix.succ.fail~w.unique, family=binomial)

fitted.yes <- aggregate(predict(fit.logit, type="response") *
                          apply(matrix.succ.fail,1,sum), by=list(W=w.cut), sum)
fitted.no <- aggregate((1-predict(fit.logit, type="response")) *
                         apply(matrix.succ.fail,1,sum), by=list(W=w.cut), sum) 
fitted.all <- matrix(c(fitted.yes$x,fitted.no$x), ncol = 2) 

#The Pearson chi-squared statistic
(x.squared = sum((observed - fitted.all)^2/fitted.all))
df <- length(observed[,1]) - length(fit.logit$coefficients)
1-pchisq(x.squared, df) 

#The likelihood ratio statistic 
glm(observed ~ seq(22.75, 29.75), family = binomial)$deviance

#Hosmer and Lemeshow Test
#install.packages("glmtoolbox")
library(glmtoolbox)
fit.grouped <- glm(observed ~ seq(22.75, 29.75), family = binomial)
hltest(fit.grouped)


#chap5c
Alcohol<-factor(c("0","<1","1-2","3-5",">=6"), levels=c("0","<1","1-2","3-5",">=6"))
malformed<-c(48,38,5,1,1)
n<-c(17066,14464,788,126,37)+malformed

#To set the first category parameter to zero, 
options(contrasts=c("contr.treatment", "contr.poly"))
(Table.5.3.logit<-glm(malformed/n~Alcohol,family=binomial, weights=n))

#More "Manually" way
# Install the required package
#install.packages("fastDummies")
# Load the library
library(fastDummies)
dta <- cbind(data.frame(Alcohol),malformed,n)
dta <- dummy_cols(data.frame(dta),select_columns = "Alcohol")
dta$Alcohol = NULL
dta$Alcohol_0 = NULL
result <- glm(malformed/n~.,family=binomial, weights=n,data=dta)

#display the result
summary(result)
cbind(logit=predict(Table.5.3.logit), fitted.prop= predict(Table.5.3.logit, type="response"))
cbind(logit=predict(result), fitted.prop= predict(result, type="response"))

#confidence interval
result.predict <- predict(result, type="response", se=T)
critval <- 1.96 ## approx 95% CI
(upr <- result.predict$fit + (critval * result.predict$se.fit))
(lwr <- result.predict$fit - (critval * result.predict$se.fit))
(fit <- result.predict$fit)


plot(result.predict$fit, axes=F, type="l",ylim = c(-0.1,0.1),
     ylab="Probability", xlab="") 
axis(2, at=seq(-0.1,0.1,.01))
lines(result.predict$fit-1.96*result.predict$se.fit,lty=3)
lines(result.predict$fit+1.96*result.predict$se.fit,lty=3)

#
(Table.5.3.logit3 <- glm(malformed/n~1, family=binomial, weights = n))


anova(Table.5.3.logit3,result,test="Chisq")

#with likelihood-ratio and Pearson chi-squared statistics
# LR statistic
summary(Table.5.3.logit3)$deviance

# Pearson chi-squared statistic
sum(residuals(Table.5.3.logit3, type="pearson")^2)

#The latter rejects the hypothesis of model fit.
1-pchisq(12.08205, df=4)


#contrast of 2-1, 3-1, 4-1, 5-1
exp(Table.5.3.logit$coefficients[2:5])
#contrast of 3-2, 4-1, 5-2
exp(Table.5.3.logit$coefficients[3:5]-Table.5.3.logit$coefficients[2])
#contrast of 4-3, 5-3
exp(Table.5.3.logit$coefficients[4:5]-Table.5.3.logit$coefficients[3])
#contrast of 5-4
exp(Table.5.3.logit$coefficients[5]-Table.5.3.logit$coefficients[4])

newdta <- c(data.frame(Alcohol = "0"),data.frame(Alcohol = "3-5"))
newdta <- data.frame(Alcohol = c("0","1-2"))
resp <- predict(Table.5.3.logit,newdta,type="response")
exp(resp[2]-resp[1])


#Linear Logit Model 
scores<-c(0,.5,1.5,4,7)
Table.5.3.LL<-glm(malformed/n~scores,family=binomial,weights=n)
summary(Table.5.3.LL)
# pearson chi-squared statistic
sum(residuals(Table.5.3.LL, type="pearson")^2)
1-pchisq(2.05229, df=3)
# LR statistic
Table.5.3.LL$null.deviance - Table.5.3.LL$deviance 
1-pchisq(4.253277, df=3)


#CAT test
x <- c(rep(scores, malformed), rep(scores, n - malformed))
y <- c(rep(1, sum(malformed)), rep(0, sum(n - malformed)))
(z2 <- 32574 * cor(x, y)^2) 
1 - pchisq(z2, df = 1) 


#chap5c
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
dta$revcolor <-  factor(dta$color, levels = c("5", "2", "3", "4")) 
dta$revspine <-  factor(dta$spine, levels = c("3", "1", "2")) 
dta$weight <- dta$weight/1000
crab.fit.logist <- glm(y ~ revcolor+ revspine+width+weight, family = binomial, data = dta)
summary(crab.fit.logist)

#Hosmer and Lemeshow Test
#install.packages("glmtoolbox")
library(glmtoolbox)
hltest(crab.fit.logist) #no lack of fit

install.packages("car")
library(car)
Anova(crab.fit.logist,type=3)

#more manually way to calculate LR chisquare
(crab.fit.logist_nocolor <- glm(y ~ revspine+width+weight, family = binomial, data = dta))
anova(crab.fit.logist_nocolor,crab.fit.logist,test="Chisq")

(crab.fit.logist_nospine <- glm(y ~ revcolor+width+weight, family = binomial, data = dta))
anova(crab.fit.logist_nospine,crab.fit.logist,test="Chisq")

(crab.fit.logist_nowidth <- glm(y ~ revspine+revcolor+weight, family = binomial, data = dta))
anova(crab.fit.logist_nowidth,crab.fit.logist,test="Chisq")



#testing without spine
Anova(crab.fit.logist_nospine,type=3)

#testing without spine or weight
crab.fit.logist_nospine_noweight <- glm(y ~ revcolor+width, family = binomial, data = dta)
Anova(crab.fit.logist_nospine_noweight,type=3)

summary(crab.fit.logist_nospine_noweight)


#Let’s see what happens when we try to combine levels of C.
#combine level 2 and 3
dta$revcolor_recode <- dta$revcolor
dta$revcolor_recode[dta$revcolor_recode == "2" | dta$revcolor_recode == "3"] = "2"
crab.fit.logist_nospine_noweight_colorrecode <- glm(y ~ revcolor_recode+width, family = binomial, data = dta)
anova(crab.fit.logist_nospine_noweight_colorrecode, crab.fit.logist_nospine_noweight, test = "Chi")
#combine all possible levels
disp.contrast <- vector()
for(i in 2:4){
  for(j in (i+1):5){
dta$revcolor_recode <- dta$revcolor
dta$revcolor_recode[dta$revcolor_recode == i | dta$revcolor_recode == j] = i
crab.fit.logist_nospine_noweight_colorrecode <- glm(y ~ revcolor_recode+width, 
                                                    family = binomial, data = dta)
result.tmp <- anova(crab.fit.logist_nospine_noweight_colorrecode, 
                    crab.fit.logist_nospine_noweight, test = "Chi")
disp.contrast <- rbind(disp.contrast,c(paste(i,"combined with",j),
                                       result.tmp$`Df`[2],result.tmp$`Deviance`[2],result.tmp$`Pr(>Chi)`[2]))
}
}

disp.contrast <- data.frame(disp.contrast)
names(disp.contrast) <- c("combination","DF","LR Chi-square","P-value")
print(disp.contrast,row.names = FALSE)


#We include dark=1; if color=4 then dark=2;
dta$darkness <- 1
dta$darkness[dta$revcolor == 5] <- 2
crab.fit.logist_nospine_noweight_darkness <- 
  glm(y ~ darkness+width, family = binomial, data = dta)
summary(crab.fit.logist_nospine_noweight_darkness)

#odds ratio
exp(coefficients(crab.fit.logist_nospine_noweight_darkness))

library(glmtoolbox)
hltest(crab.fit.logist_nospine_noweight_darkness) #no lack of fit


#Interaction model
crab.fit.logist_nospine_noweight_darkness <- 
  glm(y ~ darkness+width+darkness*width, family = binomial, data = dta)
summary(crab.fit.logist_nospine_noweight_darkness)
Anova(crab.fit.logist_nospine_noweight_darkness,type=3)
#We do not reject that the interaction is not needed.
