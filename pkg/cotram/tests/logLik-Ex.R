
## support & bounds for cotram

## libraries
library("cotram")
library("survival")

set.seed(25)

## count data
n <- 200
x <- runif(n)
y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))

## interval censored counts
yleft <- y - 1L
yleft[yleft < 0] <- -Inf
yi1 <- Surv(yleft, y, type = "interval2")


##---- log_first = FALSE ----
log_first <- FALSE

m1 <- Coxph(yi1 ~ x, bounds = b <- c(-0.01, Inf), support = s <- c(0, quantile(y, prob = .9)), log_first = log_first)
mc1 <- cotram(y ~ x, method = "cloglog", log_first = log_first)
mc1$support; mc1$bounds

L1 <- predict(m1, newdata = data.frame(yi1 = y, x = x), type = "distribution") -
  predict(m1, newdata = data.frame(yi1 = y - 1L, x = x), type = "distribution")

sum(log(L1))
logLik(m1)
logLik(mc1)

stopifnot(all.equal(log(L1), mc1$logliki(coef(as.mlt(mc1)), mc1$weights), check.attributes = FALSE))

stopifnot(all.equal(m1$logliki(coef(as.mlt(m1)), m1$weight), mc1$logliki(coef(as.mlt(mc1)), mc1$weight)))


##---- log_first = TRUE ----
log_first <- TRUE
plus_one <- as.integer(log_first)

yi2 <- Surv(yleft + plus_one, y + plus_one, type = "interval2")
m2 <- Coxph(yi2 ~ x, bounds = b + plus_one, support = s + plus_one, log_first = log_first)

mc2 <- cotram(y ~ x, method = "cloglog", log_first = log_first)
mc2$support; mc2$bounds
L2 <- predict(m2, newdata = data.frame(yi2 = y + plus_one, x = x), type = "distribution") -
  predict(m2, newdata = data.frame(yi2 = y + plus_one - 1L, x = x), type = "distribution")

sum(log(L2))
logLik(m2)
logLik(mc2)

stopifnot(all.equal(log(L2), mc2$logliki(coef(as.mlt(mc2)), mc2$weights), check.attributes = FALSE))

stopifnot(all.equal(m2$logliki(coef(as.mlt(m2)), m2$weight), mc2$logliki(coef(as.mlt(mc2)), mc2$weight)))
