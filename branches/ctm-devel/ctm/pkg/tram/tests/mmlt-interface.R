
library("tram")
library("survival")
library("sandwich")
library("numDeriv")
library("mvtnorm")

### check mmlt for wild model combinations and weights

set.seed(290875)

chk <- function(x, y, tol = 1e-3) stopifnot(all.equal(x, y, tol = tol))

N <- 100
J <- 3

### weights
w <- as.double(sample(0:5, size = N, replace = TRUE))
### X1, X2, X3 correlated
z <- matrix(rnorm(N * J), nrow = J)
L <- ltMatrices(runif(J * (J - 1) / 2))
y <- t(solve(standardize(L), z))
### uncorrelated ordered
x <- sample(gl(5, N / 5, ordered = TRUE))
### X3 right censored
e <- sample(c(TRUE, FALSE), size = N, replace = TRUE)
X3 <- exp(y[,3])
X3 <- Surv(X3, event = e)
d <- data.frame(y[,1:2], X3, x, e)
dw <- d[rep(1:nrow(d), w),,drop = FALSE]
idx <- sample(1:N, size = N / 2)

### global options (speed-up estimation)
M <- 200
args <- list(type = "ghalton", M = M)
op <- mltoptim()	### no hessian

### marginal models (use Lm to reduce number of parameters)
m1 <- Lm(X1 ~ 1, data = d)
m2 <- Lm(X2 ~ 1, data = d)
m3 <- Lm(X3 ~ 1, data = d)
m4 <- Polr(x ~ 1, data = d, method = "probit")

### joint distribution of X1, X2 (continuous), X3 (continuous and
### right-censored) and x (discrete)
m <- mmlt(m1, m2, m3, m4, data = d, args = args, optim = op)
l1 <- logLik(m)
## check gradient
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm), coef(m))
chk(s2, unname(s1))
L1 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3", "x"))

## change order of variables
m <- mmlt(m2, m4, m1, m3, data = d, args = args, optim = op)
l2 <- logLik(m)
L2 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3", "x"))

## use less evaluation points for normal integrals
m <- mmlt(m2, m4, m1, m3, data = d, optim = op,
          args = list(type = "ghalton", M = 100))
l3 <- logLik(m)
L3 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3", "x"))
## check near identity: This can't be exact, because the evaluation points
## correspond to different dimensions in the different calls
chk(l2, l1)
chk(l3, l1)
chk(L2, L1, tol = 1e-2)
chk(L3, L1, tol = 1e-2)

### continuous and ordered only
m <- mmlt(m1, m2, m4, data = d, optim = op)
l1 <- logLik(m)
l1a <- logLik(m, newdata = d[idx,]) + logLik(m, newdata = d[-idx,])
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm), coef(m))
chk(s2, unname(s1))
L1 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m2, m1, m4, data = d, optim = op)
l2 <- logLik(m)
l2a <- logLik(m, newdata = d[idx,]) + logLik(m, newdata = d[-idx,])
L2 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m4, m2, m1, data = d, optim = op)
l3 <- logLik(m)
l3a <- logLik(m, newdata = d[idx,]) + logLik(m, newdata = d[-idx,])
L3 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))
chk(l2, l1)
chk(l3, l1)
chk(L2, L1)
chk(L3, L1)

### weighted observations (1)
m1 <- Lm(X1 ~ 1, data = dw)
m2 <- Lm(X2 ~ 1, data = dw)
m3 <- Lm(X3 ~ 1, data = dw)
m4 <- Polr(x ~ 1, data = dw, method = "probit")

m <- mmlt(m1, m2, m4, data = dw, optim = op)
l1 <- logLik(m)
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm), coef(m))
chk(s2, unname(s1))
L1 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m2, m1, m4, data = dw, optim = op)
l2 <- logLik(m)
L2 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m4, m2, m1, data = dw, optim = op)
l3 <- logLik(m)
L3 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))
### check if order matters
chk(l2, l1)
chk(l3, l1)
chk(L2, L1)
chk(L3, L1)

### weighted observations (2) via weights
m1 <- Lm(X1 ~ 1, data = d, weights = w)
m2 <- Lm(X2 ~ 1, data = d, weights = w)
m3 <- Lm(X3 ~ 1, data = d, weights = w)
m4 <- Polr(x ~ 1, data = d, method = "probit", weights = w)

m <- mmlt(m1, m2, m4, data = d, optim = op)
l1w <- logLik(m)
l1wa <- logLik(m, newdata = d[idx,], w = w[idx]) + 
        logLik(m, newdata = d[-idx,], w = w[-idx])
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm, w = w), coef(m))
chk(s2, unname(s1))
L1w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m2, m1, m4, data = d, optim = op)
l2w <- logLik(m)
l2wa <- logLik(m, newdata = d[idx,], w = w[idx]) + 
        logLik(m, newdata = d[-idx,], w = w[-idx])
L2w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

m <- mmlt(m4, m2, m1, data = dw, optim = op)
l3w <- logLik(m)
l3wa <- logLik(m, newdata = d[idx,], w = w[idx]) + 
        logLik(m, newdata = d[-idx,], w = w[-idx])
L3w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "x"))

## check if order matters
chk(l2w, l1w)
chk(l3w, l1w)
chk(L2w, L1w)
chk(L3w, L1w)

## check if weights and expanded data give the same results
chk(l1w, l1wa)
chk(l2w, l2wa)
chk(l3w, l3wa)
chk(l1w, l1)
chk(l2w, l2)
chk(l3w, l3)
chk(L1w, L1)
chk(L2w, L2)
chk(L3w, L2)

### continuous and some right-censored observations
### with expanded data
m1 <- Lm(X1 ~ 1, data = dw)
m2 <- Lm(X2 ~ 1, data = dw)
m3 <- Lm(X3 ~ 1, data = dw)
m4 <- Polr(x ~ 1, data = dw, method = "probit")

m <- mmlt(m1, m2, m3, data = dw, optim = op)
l1 <- logLik(m)
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm), coef(m))
chk(s2, unname(s1))
L1 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))

m <- mmlt(m2, m1, m3, data = dw, optim = op)
l2 <- logLik(m)
L2 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))

m <- mmlt(m3, m2, m1, data = dw, optim = op)
l3 <- logLik(m)
L3 <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))
## check if order matters
chk(l2, l1)
chk(l3, l1)
chk(L2, L1)
chk(L3, L1)

### continuous and some right-censored with weights
m1 <- Lm(X1 ~ 1, data = d, weights = w)
m2 <- Lm(X2 ~ 1, data = d, weights = w)
m3 <- Lm(X3 ~ 1, data = d, weights = w)
m4 <- Polr(x ~ 1, data = d, method = "probit", weights = w)

m <- mmlt(m1, m2, m3, data = d, optim = op)
l1w <- logLik(m)
s1 <- -colSums(estfun(m))
s2 <- grad(function(parm) logLik(m, parm = parm, w = w), coef(m))
chk(s2, unname(s1))
L1w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))

m <- mmlt(m2, m1, m3, data = d, optim = op)
l2w <- logLik(m)
L2w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))

m <- mmlt(m3, m2, m1, data = d, optim = op)
l3w <- logLik(m)
L3w <- aperm(coef(m, type = "Lambda"), perm = c("X1", "X2", "X3"))

## check if order matters
chk(l2w, l1w)
chk(l3w, l1w)
chk(L2w, L1w)
chk(L3w, L1w)

## check if weights and expanded data give the same results
chk(l1w, l1)
chk(l2w, l2)
chk(l3w, l3)
chk(L1w, L1)
chk(L2w, L2)
chk(L3w, L2)

### simulation from fitted / unfitted models
data("iris")

m1 <- BoxCox(Sepal.Length ~ Species, data = iris)
m2 <- BoxCox(Petal.Length ~ Species, data = iris)

### fit model
m <- mmlt(m1, m2, data = iris)
cf <- coef(m)

### resample data
s1 <- as.data.frame(simulate(m, newdata = iris, seed = 290875))

### change coeffcient (to zero correlation)
cf[length(cf)] <- 0
coef(m) <- cf

s2 <- as.data.frame(simulate(m, newdata = iris, seed = 290875))

### expect _different_ results
stopifnot(!isTRUE(all.equal(s1, s2)))

### set-up model shell
m0 <- mmlt(m1, m2, data = iris, theta = cf, dofit = FALSE)

s3 <- as.data.frame(simulate(m, newdata = iris, seed = 290875))
### these are identical
chk(s2, s3)
