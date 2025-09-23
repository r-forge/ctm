
library("basefun")

### check first and second derivatives
if (require("numDeriv")) {

OR <- 5
cf <- cumsum(runif(OR + 1))
### exclude support b/c of lack of differentiability
(y <- c(1:9, 11:19, 21:29) / 10)

test <- expression({
  yv <- numeric_var("y", support = c(1, 2), bounds = c(.1, 2.9))
  bs1 <- Bernstein_basis(yv, log_first = LF, order = OR, extrapolate = EX)
  h1 <- model.matrix(bs1, data = data.frame(y = y)) %*% cf
  h1p <- model.matrix(bs1, data = data.frame(y = y), deriv = c("y" = 1)) %*% cf
  h1pp <- model.matrix(bs1, data = data.frame(y = y), deriv = c("y" = 2)) %*% cf

  fun <- function(x)
    model.matrix(bs1, data = data.frame(y = x), deriv = c("y" = 0)) %*% cf

  stopifnot(all.equal(sapply(y, function(x) grad(fun, x)), c(h1p),
                      tol = 1e-5))

  fun <- function(x)
    model.matrix(bs1, data = data.frame(y = x), deriv = c("y" = 1)) %*% cf

  stopifnot(all.equal(sapply(y, function(x) grad(fun, x)), c(h1pp),
                      tol = 1e-5))
})

LF <- FALSE
EX <- FALSE

eval(test)

LF <- TRUE
EX <- FALSE

eval(test)

LF <- FALSE
EX <- TRUE

eval(test)

LF <- TRUE
EX <- TRUE

eval(test)

}

