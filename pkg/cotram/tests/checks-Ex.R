## cotram: test checks and settings

library("cotram")
set.seed(25)

y <- 0L:100L
x <- runif(length(y))

yp1 <- y + 1L
d <- data.frame(y = y, x = x)

m <- cotram(y ~ x, log_first = TRUE)

mm <- cotram(y ~ x, log_first = FALSE)

## ---- checks for y > 0 and y %% 1 == 0 ----
ym <- -y
try(cotram(ym ~ x))
yn <- y/2
try(cotram(yn ~ x))

try(logLik(mm, newdata = data.frame(y = ym, x = x)))
try(logLik(mm, newdata = data.frame(y = yn, x = x)))

try(predict(mm, newdata = data.frame(y = ym, x = x)))
try(predict(mm, newdata = data.frame(y = yn, x = x)))

try(plot(mm, newdata = data.frame(y = ym, x = x)))
try(plot(mm, newdata = data.frame(y = yn, x = x)))

try(plot(mm, q = ym, newdata = model.frame(mm)))
try(plot(mm, q = yn, newdata = model.frame(mm)))
