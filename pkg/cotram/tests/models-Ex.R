## comparing log-likelihoods and coefficients

library("cotram")
library("survival")

set.seed(29)

## dgp
n <- 200
x <- runif(n)
y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))
yn <- as.numeric(y)

## interval censored counts
yleft <- y - 1L
yleft[yleft < 0] <- -Inf
yi1 <- Surv(yleft, y, type = "interval2")
yip1 <- Surv(yleft + 1L, y + 1L, type = "interval2")

df <- data.frame(x = x, y = y, yn = yn, yip1 = yip1)

trainID <- sample(1:nrow(df), size = 0.25 * nrow(df))
d <- df[trainID,]
nd <- df[-trainID,]

## test model
m1 <- cotram(y ~ x, data = d, method = "cloglog", 
              log_first = TRUE)

m1n <- cotram(yn ~ x, data = d, method = "cloglog",
              log_first = TRUE)

m2 <- Coxph(yip1 ~ x, data = d, support = m1$support, bounds = m1$bounds,
            log_first = TRUE)

m3 <- cotram(y ~ x , data = d, method = "cloglog",
              log_first = FALSE)


## compare coefficients
cf1 <- coef(as.mlt(m1))
# cf1n <- coef(as.mlt(m1n))
# stopifnot(all.equal(cf1n, cf1, check.attributes = FALSE))

c2 <- as.mlt(m2)
cf2 <- coef(c2)
cf2[c2$shiftcoef] <- -cf2[c2$shiftcoef]
stopifnot(all.equal(cf1, cf2, check.attributes = FALSE))

# compare likelihoods
l1 <- m1$logliki(coef(as.mlt(m1)), rep(1, length(d$y)))
l1n <- m1$logliki(coef(as.mlt(m1n)), rep(1, length(d$yn)))
stopifnot(all.equal(l1, l1n))

l2 <- m2$logliki(coef(as.mlt(m2)), rep(1, length(d$yp1)))
stopifnot(all.equal(l1, l2))

# logLik for newdata
stopifnot(isTRUE(logLik(m1) == logLik(m1, newdata = d)))
stopifnot(isTRUE(logLik(m1, newdata = nd) == logLik(m2, newdata = nd)))
stopifnot(isTRUE(logLik(m1, newdata = nd) == logLik(m1n, newdata = nd)))
