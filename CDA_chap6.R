#chap6a
dta <- read.table("Horseshoe crab data.txt", header = TRUE, sep = "")
dta$revcolor <-  factor(dta$color, levels = c("5", "2", "3", "4")) 
dta$revspine <-  factor(dta$spine, levels = c("3", "1", "2")) 
dta$weight <- dta$weight/1000
# LR statistic
crab.fit.logist$null.deviance - crab.fit.logist$deviance 
1-pchisq(40.55652, df=7)


# create a large pool of potential predictors
predictors <- c("revcolor", "revspine", "width", "weight", 
                "revcolor:revspine", "revcolor:width", "revcolor:weight",
                "revspine:width", "revspine:weight",
                "width:weight")
formula_str <- paste("y ~", paste(predictors, collapse = "+"))
model_full <- glm(formula_str, data = dta, family = binomial)
# more manually perform backward elimination with cutoff based on Wald p-values
p_cutoff <- 0.15
while (TRUE) {
#---------- This is for LRT -----------#
  result_full <- anova(model_full, test = "LRT")
  p_values <- result_full$`Pr(>Chi)`
  p_values <- p_values[-1] #remove the NA for the intercept
#---------- This is for Wald ---------#
#  result_full <- Anova(model_full,test.statistic = "Wald") #using "car" package
#  p_values <- result_full$`Pr(>Chisq)`

  
  
  # identify predictors with p-values larger than the cutoff
  remove_predictors <- predictors[p_values > p_cutoff]
  
  # stop the loop if no more predictors to remove
  if (length(remove_predictors) == 0) {
    break
  }else{
    remove_predictors <- predictors[p_values == max(p_values[is.na(p_values)==FALSE])]
  }
 # print(paste(remove_predictors,result_full$Df[p_values == max(p_values[is.na(p_values)==FALSE])],max(p_values[is.na(p_values)==FALSE]),sep=" "))
  # remove the predictors with large p-values
  print(remove_predictors)
  predictors <- setdiff(predictors, remove_predictors)
  
  # refit the model with the reduced set of predictors
  formula_str <- paste("y ~", paste(predictors, collapse = "+"))
#  print(formula_str)
  model_full <- glm(formula_str, data = dta, family = binomial)
}
print(summary(model_full))
