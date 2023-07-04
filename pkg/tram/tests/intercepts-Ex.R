
library("tram")
options(digits = 1) ## make M1 happy

## dgp
set.seed(29)
N <- 1000
x <- runif(N)
y <- rchisq(N, df = 4)
d <- data.frame(y = y, x = x, one = 1)

d$xx <- d$x - 5 ## shift x


## no h-intercept; additional mu-intercept
m0 <- Survreg(y ~ x | x, data = d, remove_intercept = FALSE)
coef(as.mlt(m0))
sqrt(diag(vcov(m0)))

m1 <- Survreg(y ~ xx | x, data = d, remove_intercept = FALSE)
coef(as.mlt(m1))
sqrt(diag(vcov(m1)))

m2 <- Survreg(y ~ xx | xx, data = d, remove_intercept = FALSE)
coef(as.mlt(m2))
sqrt(diag(vcov(m2)))


## h-intercept; scale_shift
m0 <- Survreg(y ~ x | x, data = d, scale_shift = TRUE)
coef(as.mlt(m0))
sqrt(diag(vcov(m0)))

m1 <- Survreg(y ~ xx | x, data = d, scale_shift = TRUE)
coef(as.mlt(m1))
sqrt(diag(vcov(m1)))

m2 <- Survreg(y ~ xx | xx, data = d, scale_shift = TRUE)
coef(as.mlt(m2))
sqrt(diag(vcov(m2)))

OR <- 1

d$ly <- log(y)

## no h-intercept; additional mu-intercept
m0 <- Coxph(ly ~ x | x, data = d, log_first = FALSE, order = OR, remove_intercept = FALSE)
coef(as.mlt(m0))
sqrt(diag(vcov(m0)))

m1 <- Coxph(ly ~ xx | x, data = d, log_first = FALSE, order = OR, remove_intercept = FALSE)
coef(as.mlt(m1))
sqrt(diag(vcov(m1)))

m2 <- Coxph(ly ~ xx | xx, data = d, log_first = FALSE, order = OR, remove_intercept = FALSE)
coef(as.mlt(m2))
sqrt(diag(vcov(m2)))

OR <- 1

## h-intercept; scale_shift
m0 <- Coxph(y ~ x | x, data = d, scale_shift = TRUE, log_first = TRUE, order = OR)
coef(as.mlt(m0))
sqrt(diag(vcov(m0)))

m1 <- Coxph(y ~ xx | x, data = d, scale_shift = TRUE, log_first = TRUE, order = OR)
coef(as.mlt(m1))
sqrt(diag(vcov(m1)))

m2 <- Coxph(y ~ xx | xx, data = d, scale_shift = TRUE, log_first = TRUE, order = OR)
coef(as.mlt(m2))
sqrt(diag(vcov(m2)))

