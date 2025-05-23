
library("basefun")
n <- 100
set.seed(29)

### one-dim problem
## set-up simple regression problem
x <- sort(runif(100))
y <- x^2 + rnorm(length(x), sd = .1) 
## set-up monotone Bernstein polynom basis
xvar <- numeric_var("x", support = c(0.0, 1.0))
Bb <- Bernstein_basis(xvar, order = 5, ui = "increasing")
## evaluate basis
X1 <- model.matrix(Bb, data = data.frame(x = x))
X2 <- Bb(x)
stopifnot(isTRUE(all.equal(X1, X2)))
## fit model
m1 <- lm(y ~  X1 - 1)
m2 <- lm(y ~  Bb(x) - 1, data = data.frame(y = y, x = x))
Bb0 <- Bernstein_basis(xvar, order = 5, ui = c("increasing", "zeroint"))
X0 <- model.matrix(Bb0, data = data.frame(x = x))
m0 <- lm(y ~  X0)
stopifnot(isTRUE(all.equal(fitted(m1), fitted(m2))))
stopifnot(isTRUE(all.equal(fitted(m1), fitted(m0))))
stopifnot(isTRUE(all.equal(sum(coef(m1)) / length(coef(m1)), coef(m0)["(Intercept)"], 
                    check.attributes = FALSE)))

if (require("quadprog")) {
  ### check contraints fits
  q1 <- solve.QP(crossprod(X1), crossprod(X1, y), 
                 t(as(attr(X1, "constraint")$ui, "matrix")), attr(X1, "constraint")$ci)
  q0 <- solve.QP(crossprod(cbind(1, X0)), crossprod(cbind(1, X0), y), 
                 t(cbind(0, as(attr(X0, "constraint")$ui, "matrix"))), attr(X0, "constraint")$ci)
  stopifnot(isTRUE(all.equal(X1 %*% q1$solution, cbind(1, X0) %*% q0$solution)))
  stopifnot(isTRUE(all.equal(sum(q1$solution) / length(q1$solution), q0$solution[1L])))
  ### check derivatives
  stopifnot(isTRUE(all.equal(model.matrix(Bb, data = data.frame(x = x), deriv = c("x" = 1)) %*% q1$solution,
                             model.matrix(Bb0, data = data.frame(x = x), deriv = c("x" = 1)) %*% q0$solution[-1L])))
}



stopifnot(isTRUE(all.equal(coef(m1), coef(m2), check.attributes = FALSE)))
## generate new data from support of basis
xn <- mkgrid(Bb, n = 100)
## compute estimated regression function
p1 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1))
p2 <- predict(m2, newdata = data.frame(x = xn)) 
stopifnot(isTRUE(all.equal(c(p1), p2, check.attributes = FALSE)))
## compute derivative of estimated regression function
dp1 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1), deriv = c(x = 1))
dp0 <- predict(Bb0, newdata = data.frame(x = xn), coef = coef(m0)[-1], deriv = c(x = 1))
stopifnot(isTRUE(all.equal(dp1, dp0)))
dp12 <- predict(Bb, newdata = data.frame(x = xn), coef = coef(m1), deriv = c(x = 1), integrate = TRUE)
unique(round(c(p1 - dp12), 5)) ### the same up to a constant

### two-dim (ANCOVA) problem
gf <- gl(3, 1)
g <- sample(gf, length(x), replace = TRUE)
y <- x^2 + (g == "2") * sin(x) + (g == "3") * cos(x) + rnorm(length(x), sd = .05)
## generate a basis for a factor (treatment contrasts)
## this is equal to model.matrix(~ gf)
gb <- as.basis(~ g, remove_intercept = FALSE, data = data.frame(g = gf))
## join the two bases by the kronecker product
bb <- b(b1 = Bb, b2 = gb)
## evaluate new two-dim basis
X1 <- model.matrix(bb, data = data.frame(x = x, g = g))
## fit model
m1 <- lm(y ~  X1 - 1)
m2 <- lm(y ~  model.matrix(bb, data.frame(x = x, g = g)) - 1, data = data.frame(y = y, x = x, g = g))
stopifnot(isTRUE(all.equal(coef(m1), coef(m2), check.attributes = FALSE)))
## compute estimated regression functions
d <- mkgrid(bb, n = 100)
## for each group
p1 <- sapply(gf, function(l) predict(bb, newdata = data.frame(x = d$x, g = l), coef = coef(m1)))
## the same via _linear array_ approach
p2 <- predict(bb, newdata = d, coef(m1))
## brute force; 2 times
p3 <- matrix(predict(bb, newdata = do.call(expand.grid, d), coef(m1)), ncol = nlevels(gf))
p4 <- matrix(predict(m2, newdata = do.call(expand.grid, d)), ncol = nlevels(gf))
stopifnot(isTRUE(all.equal(p1, p2, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(p2, p3, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(p3, p4, check.attributes = FALSE)))
## compute derivative wrt the first element
dp2 <- predict(bb, newdata = d, coef(m1), deriv = c(x = 1))


## join the two bases additively
gb <- as.basis(~ g, remove_intercept = TRUE, data = data.frame(g = gf))
bb <- c(b1 = Bb, b2 = gb)
## evaluate new two-dim basis
X1 <- model.matrix(bb, data = data.frame(x = x, g = g))
## fit model
m1 <- lm(y ~  X1 - 1)
m2 <- lm(y ~  model.matrix(bb, data.frame(x = x, g = g)) - 1, data = data.frame(y = y, x = x, g = g))
stopifnot(isTRUE(all.equal(coef(m1), coef(m2), check.attributes = FALSE)))
## compute estimated regression functions
d <- mkgrid(bb, n = 100)
## for each group
p1 <- c(sapply(gf, function(l) predict(bb, newdata = data.frame(x = d$x, g = l), coef = coef(m1))))
## the same via expand.grid approach
p2 <- c(predict(bb, newdata = d, coef(m1)))
## brute force; 2 times
p3 <- predict(bb, newdata = do.call(expand.grid, d), coef(m1))
p4 <- predict(m2, newdata = do.call(expand.grid, d))
stopifnot(isTRUE(all.equal(p1, p2, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(p2, p3, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(p3, p4, check.attributes = FALSE)))
## compute derivative wrt the first element
dp2 <- predict(bb, newdata = d, coef(m1), deriv = c(x = 1))

### a third variable
z <- runif(length(x))
zvar <- numeric_var("z", support = c(0.0, 1.0))
bz <- as.basis(~ z - 1, data = zvar)

testb <- function(bb) {
    X1 <- model.matrix(bb, data = data.frame(x = x, g = g, z = z))
    ## fit model
    m1 <- lm(y ~  X1 - 1)
    m2 <- lm(y ~  model.matrix(bb, data.frame(x = x, g = g, z = z)) - 1, 
             data = data.frame(y = y, x = x, g = g, z = z))
    stopifnot(isTRUE(all.equal(coef(m1), coef(m2), check.attributes = FALSE)))
    ## compute estimated regression functions
    d <- mkgrid(bb, n = 100)
    p2 <- c(predict(bb, newdata = d, coef(m1)))
    ## brute force; 2 times
    p3 <- predict(bb, newdata = do.call(expand.grid, d), coef(m1))
    p4 <- predict(m2, newdata = do.call(expand.grid, d))
    stopifnot(isTRUE(all.equal(p2, p3, check.attributes = FALSE)))
    stopifnot(isTRUE(all.equal(p3, p4, check.attributes = FALSE)))
    ## compute derivative wrt the first element
    dp2 <- predict(bb, newdata = d, coef(m1), deriv = c(x = 1))
}

testb(c(b1 = Bb, b2 = b(b2g = gb, b2z = bz)))
testb(c(b1 = Bb, b2 = c(b2g = gb, b2z = bz)))
testb(b(b1 = Bb, b2 = c(b2g = gb, b2z = bz)))

testb(c(b0 = c(b1 = Bb, b2 = b(b2g = gb, b2z = bz))))
testb(c(b0 = c(b1 = Bb, b2 = c(b2g = gb, b2z = bz))))
testb(c(b0 = b(b1 = Bb, b2 = c(b2g = gb, b2z = bz))))

### cyclic basis
cb <- cyclic_basis(numeric_var("x"), order = 3, frequency = pi)
### generate data + coefficients
x <- data.frame(x = -10:10 * pi)
f <- cb(x) %*% runif(6)
stopifnot(all(abs(diff(f)) < sqrt(.Machine$double.eps)))

