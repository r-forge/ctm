
library("cotram")
options(digits = 2)

set.seed(25)

## dgp
n <- 200
x <- runif(n, 0, 20)
y <- as.integer(rnbinom(n, mu = exp(.2 + .1 * x), size = 3))
d <- data.frame(y = y, x = x)

m <- cotram(y ~ x, data = d)

nd <- model.frame(m)[3,]

## Confband for grid of counts
confband(m, type = "distribution", newdata = nd)

## Confband for K grid points
confband(m, type = "distribution", newdata = nd, smooth = TRUE, K = 40)


if (FALSE){
layout(matrix(1:2, nrow = 1))
type = "trafo"

cb <- confband(m, type =  type, newdata = nd)
plot(m, type = type, newdata = nd, 
     confidence = "band", col = "red", ylim = c(-2, 15))
lines(x = cb[, "q"], y = cb[, "lwr"], type = "s")
lines(x = cb[, "q"], y = cb[, "upr"], type = "s")

cb.s <- confband(m, type = type, newdata = nd, smooth = TRUE)
plot(m, type = type, newdata = nd, smooth = TRUE,
     confidence = "band", col = "red", ylim = c(-2, 15))
lines(x = cb.s[, "q"], y = cb.s[, "lwr"], type = "l")
lines(x = cb.s[, "q"], y = cb.s[, "upr"], type = "l")
}

