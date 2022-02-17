## test checks and settings

library("cotram")
set.seed(25)

y <- 0L:100L
x <- runif(length(y))

yp1 <- y + 1L
d <- data.frame(y = y, x = x)

m <- cotram(y ~ x, log_first = TRUE)

mm <- cotram(y ~ x, log_first = FALSE)

## ---- checks for y > 0 and y %% 1 == 0 ----

## quick check that it returns a error for non-positive & non-integers
.check_error <- function(expr) stopifnot(class(try(expr)) == "try-error")

## negative y
ym <- -y
.check_error(cotram(ym ~ x))

## non-integer y
yn <- y/2
.check_error(cotram(yn ~ x))

.check_error(logLik(mm, newdata = data.frame(y = ym, x = x)))
.check_error(logLik(mm, newdata = data.frame(y = yn, x = x)))

.check_error(predict(mm, newdata = data.frame(y = ym, x = x)))
.check_error(predict(mm, newdata = data.frame(y = yn, x = x)))

.check_error(plot(mm, newdata = data.frame(y = ym, x = x)))
.check_error(plot(mm, newdata = data.frame(y = yn, x = x)))

.check_error(confband(mm, q = ym, newdata = model.frame(mm)))
.check_error(confband(mm, q = yn, newdata = model.frame(mm)))

.check_error(plot(mm, q = ym, newdata = model.frame(mm)))
.check_error(plot(mm, q = yn, newdata = model.frame(mm)))
