set.seed(123)  # for reproducibility

# Generate two normally distributed covariates
n <- 1000
x1 <- rnorm(n, 0, 1)
x2 <- rnorm(n, 0, 1)

# Generate binary response variable with complete separation
y <- ifelse(x1 + x2 > 0, 1, 0)

# Fit logistic regression model separately
fit <- glm(y ~ x1, family = binomial(link = "logit"))
summary(fit)
fit <- glm(y ~ x2, family = binomial(link = "logit"))
summary(fit)

# Fit logistic regression model
fit <- glm(y ~ x1 + x2, family = binomial(link = "logit"))
summary(fit)

plot(x1, y)
plot(x2, y)
plot(x1+x2, y)
