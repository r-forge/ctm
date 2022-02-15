## cotram: test checks and settings

library("cotram")
set.seed(25)

y <- 0L:100L
x <- runif(length(y))

yp1 <- y + 1L
mt <- Coxph(yp1 ~ x)

m <- cotram(y ~ x)

mm <- cotram(y ~ x, log_first = FALSE)

nd <- data.frame(y = y, x = x)

## ---- checks for y > 0 and y %% 1 == 0 ----
ym <- -y
try(cotram(ym ~ x))
yn <- y/2
try(cotram(yn ~ x))

try(logLik(m, newdata = data.frame(y = ym, x = x)))
try(logLik(m, newdata = data.frame(y = yn, x = x)))

try(predict(m, newdata = data.frame(y = ym, x = x)))
try(predict(m, newdata = data.frame(y = yn, x = x)))

try(plot(m, newdata = data.frame(y = ym, x = x)))
try(plot(m, newdata = data.frame(y = yn, x = x)))

try(plot(m, q = ym, newdata = model.frame(m, x = x)))
try(plot(m, q = yn, newdata = model.frame(m, x = x)))

## --- checks for plus_one ----
stopifnot(all(model.frame(m)["y"] == y)) ## model.frame of y for log_first

stopifnot(all(mkgrid(m)$y == y))

stopifnot(logLik(m) == logLik(m, newdata = nd))
stopifnot(logLik(mm) == logLik(mm, newdata = nd))

## predict
ppr <- predict(mm, newdata = nd, type = "density")
pr <- predict(m, newdata = nd, type = "density")

stopifnot(rownames(pr) == rownames(ppr))
stopifnot(!any(abs(ppr - pr) > 0.1))

pprq <- predict(mm, newdata = nd, type = "density", q = 0:10)
prq <- predict(m, newdata = nd, type = "density", q = 0:10)

stopifnot(rownames(prq) == rownames(pprq))
stopifnot(!any(abs(pprq - prq) > 0.1))

## plot
if (FALSE) {
plot(m, newdata = nd,  type = "density", col = 1)
plot(mm, newdata = nd, type = "density", col = 2, add = TRUE)
abline(v = 0)

plot(m, newdata = nd,  type = "logdensity", col = 1, ylim = c(-6, 0))
plot(mm, newdata = nd, type = "logdensity", col = 2, add = TRUE, ylim = c(-6, 0))
abline(v = 0)

plot(m, newdata = nd, type = "distribution", col = 1)
plot(mm, newdata = nd, type = "distribution", col = 2, add = TRUE)
abline(v = 0)

plot(m, newdata = nd, type = "trafo", col = 1)
plot(mm, newdata = nd, type = "trafo", col = 2, add = TRUE)   
abline(v = 0)

plot(m, newdata = nd, type = "survivor", col = 1)
plot(mm, newdata = nd, type = "survivor", col = 2, add = TRUE)   
abline(v = 0)

plot(m, newdata = nd, type = "cumhazard", col = 1)
plot(mm, newdata = nd, type = "cumhazard", col = 2, add = TRUE)   
abline(v = 0)
}