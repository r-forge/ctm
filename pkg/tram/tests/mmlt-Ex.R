
library("tram")
library("mvtnorm")

set.seed(25)
chk <- function(...) all.equal(...)

J <- 4
N <- 1000
S <- cov2cor(tcrossprod(matrix(runif(J * J), ncol = J)))
y <- rmvnorm(N, sigma = S)
u <- as.data.frame(plogis(y))
x <- runif(N)
d <- cbind(u, x)
un <- colnames(d)[1:J]

m <- lapply(un, function(i)
    BoxCox(as.formula(paste(i, "~ x")), data = d, bounds = c(0, 1), support = c(0, 1)))
m$data <- d
m$formula <- ~ 1
mm <- do.call("mmlt", m)

chk(c(logLik(mm)), sum(predict(mm, newdata = d, type = "density", log = TRUE)))
L <- as.array(coef(mm, type = "Lambda"))[,,1]
chk(as.array(coef(mm, type = "Lambdainv"))[,,1], solve(L))
chk(as.array(coef(mm, type = "Sigma"))[,,1], tcrossprod(solve(L)))
chk(as.array(coef(mm, type = "Cor"))[,,1], cov2cor(tcrossprod(solve(L))))

### marginal normal
m$conditional <- FALSE
mmN <- do.call("mmlt", m)

chk(logLik(mm), logLik(mmN))
chk(c(logLik(mmN)), sum(predict(mmN, newdata = d, type = "density", log = TRUE)))

chk(as.array(coef(mm, type = "Lambda"))[,,1], 
    as.array(coef(mmN, type = "Lambda"))[,,1])
chk(as.array(coef(mm, type = "Lambdainv"))[,,1], 
    as.array(coef(mmN, type = "Lambdainv"))[,,1])
chk(as.array(coef(mm, type = "Sigma"))[,,1], 
    as.array(coef(mmN, type = "Sigma"))[,,1])
chk(as.array(coef(mm, type = "Spearman"))[,,1], 
    as.array(coef(mmN, type = "Spearman"))[,,1])

chk(predict(mm, newdata = d, type = "density", log = TRUE), 
    predict(mmN, newdata = d, type = "density", log = TRUE))
chk(predict(mm, newdata = d, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, type = "distribution", log = TRUE))

chk(predict(mm, newdata = d, margins = 1:2, type = "density", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:2, type = "density", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:2, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:2, type = "distribution", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:3, type = "density", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:3, type = "density", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:3, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:3, type = "distribution", log = TRUE))

chk(sapply(1:J, function(i) predict(mm, margins = i, newdata = d, type = "density", log = TRUE)),
    sapply(1:J, function(i) predict(mmN, margins = i, newdata = d, type = "density", log = TRUE)), 
    check.attributes = FALSE)

m <- lapply(un, function(i)
    Colr(as.formula(paste(i, "~ x")), data = d, bounds = c(0, 1), support = c(0, 1)))
m$data <- d
m$formula <- ~ 1
mmC <- do.call("mmlt", m)

chk(c(logLik(mmC)), sum(predict(mmC, newdata = d, type = "density", log = TRUE)))
logLik(mmC)

m <- lapply(un, function(i)
    BoxCox(as.formula(paste(i, "~ x")), data = d, bounds = c(0, 1), support = c(0, 1)))
m$data <- d
m$formula <- ~ x
mm <- do.call("mmlt", m)

chk(c(logLik(mm)), sum(predict(mm, newdata = d, type = "density", log = TRUE)))
L <- as.array(coef(mm, newdata = d, type = "Lambda"))[,,1]
chk(as.array(coef(mm, newdata = d, type = "Lambdainv"))[,,1], solve(L))
chk(as.array(coef(mm, newdata = d, type = "Sigma"))[,,1], tcrossprod(solve(L)))
chk(as.array(coef(mm, newdata = d, type = "Cor"))[,,1], cov2cor(tcrossprod(solve(L))))

### fake normal
for (j in 1:J) m[[j]]$todistr$name <- "CarlFriedrich"

mmN <- do.call("mmlt", m)

chk(logLik(mm), logLik(mmN))
chk(c(logLik(mmN)), sum(predict(mmN, newdata = d, type = "density", log = TRUE)))

chk(as.array(coef(mm, newdata = d, type = "Lambda"))[,,1], 
    as.array(coef(mmN, newdata = d, type = "Lambda"))[,,1])
chk(as.array(coef(mm, newdata = d, type = "Lambdainv"))[,,1], 
    as.array(coef(mmN, newdata = d, type = "Lambdainv"))[,,1])
chk(as.array(coef(mm, newdata = d, type = "Sigma"))[,,1], 
    as.array(coef(mmN, newdata = d, type = "Sigma"))[,,1])
chk(as.array(coef(mm, newdata = d, type = "Spearman"))[,,1], 
    as.array(coef(mmN, newdata = d, type = "Spearman"))[,,1])

chk(predict(mm, newdata = d, type = "density", log = TRUE), 
    predict(mmN, newdata = d, type = "density", log = TRUE))
chk(predict(mm, newdata = d, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, type = "distribution", log = TRUE))

chk(predict(mm, newdata = d, margins = 1:2, type = "density", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:2, type = "density", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:2, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:2, type = "distribution", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:3, type = "density", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:3, type = "density", log = TRUE))
chk(predict(mm, newdata = d, margins = 1:3, type = "distribution", log = TRUE), 
    predict(mmN, newdata = d, margins = 1:3, type = "distribution", log = TRUE))

chk(sapply(1:J, function(i) predict(mm, margins = i, newdata = d, type = "density", log = TRUE)),
    sapply(1:J, function(i) predict(mmN, margins = i, newdata = d, type = "density", log = TRUE)), 
    check.attributes = FALSE)

m <- lapply(un, function(i)
    Colr(as.formula(paste(i, "~ x")), data = d, bounds = c(0, 1), support = c(0, 1)))
m$data <- d
m$formula <- ~ x
mmC <- do.call("mmlt", m)

chk(c(logLik(mmC)), sum(predict(mmC, newdata = d, type = "density", log = TRUE)))
logLik(mmC)

##### FIRST SCENARIO: CONSTANT LAMBDA #####
set.seed(29)

ll <- numeric(50)

N <- 5000
p <- 3
X <- matrix(runif(N * p), ncol = p)
m1 <- 1 + X %*% c(2, 1, 1)
m2 <- 1 + X %*% c(1, 2, 1)
lb <- (off <- 0.5) + X %*% (cf <- c(0, 2, 0))
d <- data.frame(X)
Y <- matrix(NA, nrow = N, ncol = 2)
colnames(Y) <- c("Y1", "Y2")

cr <- numeric(N)
for (i in 1:N) {
  Si <- diag(2)
  Si[1,2] <- Si[2,1] <- .5
  cr[i] <- cov2cor(Si)[2,1]
  
  Y[i,] <- rmvnorm(1, mean = c(m1[i], m2[i]), sigma = Si)
}


##### only BoxCox margins: ##### 
d <- cbind(d, Y)
b1 <- as.mlt(Lm(Y1 ~ X1 + X2 + X3, data = d))
b2 <- as.mlt(Lm(Y2 ~ X1 + X2 + X3, data = d))

## constant correlations. expect identical logliks and lambda parameters
mm01 <- mmlt(b1, b2, formula = ~ 1, data = d)
mm02 <- mmlt(b2, b1, formula = ~ 1, data = d)

logLik(mm01) 
logLik(mm02)

coef(mm01)["Y2.Y1.(Intercept)"]
coef(mm02)["Y1.Y2.(Intercept)"]

## checking gradients
all.equal(c(numDeriv::grad(mm01$ll, mm02$par)),c(mm01$sc(mm02$par)), 
          check.attributes = FALSE, tol = 1e-4)	

## predicting marginal distributions and comparing across models with constant lambda
predict(mm01, newdata = d[1:5,], q = -2:2, 
        margins = 1, type = "distribution")
predict(mm02, newdata = d[1:5,], q = -2:2, 
        margins = 2, type = "distribution")

## expect correlations to be the same for the model with constant lambdas
c(coef(mm01, newdata = d[1:5,], type = "Cor"))
c(coef(mm02, newdata = d[1:5,], type = "Cor"))


##### mix of BoxCox and Colr margins: ##### 
d$Y1 <- (d$Y1 - min(d$Y1))/(max(d$Y1) - min(d$Y1))

b1 <- as.mlt(Colr(Y1 ~ X1 + X2 + X3, data = d, order = 1))
b2 <- as.mlt(Lm(Y2 ~ X1 + X2 + X3, data = d))

mm01 <- mmlt(b1, b2, formula = ~ 1, data = d)
mm02 <- mmlt(b2, b1, formula = ~ 1, data = d)

logLik(mm01)
logLik(mm02)

coef(b1)
coef(b2)
coef(mm01)
coef(mm02)

## checking gradient
all.equal(c(numDeriv::grad(mm01$ll, coef(mm01))), c(mm01$sc(coef(mm01))), 
          check.attributes = FALSE, tol = 1e-4)
all.equal(c(numDeriv::grad(mm02$ll, coef(mm02))), c(mm02$sc(coef(mm02))), 
          check.attributes = FALSE, tol = 1e-4)


##### SECOND SCENARIO: COVARIATE DEPENDENT LAMBDA #####
set.seed(29)

ll <- numeric(50)

N <- 5000
p <- 3
X <- matrix(runif(N * p), ncol = p)
m1 <- 1 + X %*% c(2, 1, 1)
m2 <- 1 + X %*% c(1, 2, 1)
lb <- (off <- 0.5) + X %*% (cf <- c(0, 2, 0))
d <- data.frame(X)
Y <- matrix(NA, nrow = N, ncol = 2)
colnames(Y) <- c("Y1", "Y2")

cr <- numeric(N)
for (i in 1:N) {
  L <- diag(2)
  L[2,1] <- lb[i]
  Si <- solve(L) %*% t(solve(L))
  cr[i] <- cov2cor(Si)[2,1]
  
  Y[i,] <- rmvnorm(1, mean = c(m1[i], m2[i]), sigma = Si)
}


##### only BoxCox margins: ##### 
d <- cbind(d, Y)
b1 <- as.mlt(Lm(Y1 ~ X1 + X2 + X3, data = d))
b2 <- as.mlt(Lm(Y2 ~ X1 + X2 + X3, data = d))

## constant correlations. expect identical logliks and lambda parameters
mm01 <- mmlt(b1, b2, formula = ~ 1, data = d)
mm02 <- mmlt(b2, b1, formula = ~ 1, data = d)

logLik(mm01) 
logLik(mm02)

## x-dependent correlations. expect slightly different logliks
mm1 <- mmlt(b1, b2, formula = ~ X1 + X2 + X3, data = d)
mm2 <- mmlt(b2, b1, formula = ~ X1 + X2 + X3, data = d)

logLik(mm1)
logLik(mm2)


## checking gradients
all.equal(c(numDeriv::grad(mm01$ll, mm02$par)),c(mm01$sc(mm02$par)), 
          check.attributes = FALSE, tol = 1e-4)
all.equal(c(numDeriv::grad(mm1$ll, mm2$par)),c(mm1$sc(mm2$par)),
          check.attributes = FALSE, tol = 1e-4)

## plotting correlations
x <- 0:4 / 4
nd <- expand.grid(X1 = x, X2 = x, X3 = x)
lhat <- off + as.matrix(nd) %*% cf
chat <- sapply(lhat, function(x) {
  L <- diag(2)
  L[2,1] <- x
  Si <- solve(L) %*% t(solve(L))
  cov2cor(Si)[2,1]
})
CR <- cbind(chat, coef(mm1, newdata = nd, type = "Cor"), 
            coef(mm2, newdata = nd, type = "Cor"))
pairs(CR)

## predicting marginal distributions and comparing across models with constant lambda
predict(mm01, newdata = nd[1:5,], q = -2:2, 
        margins = 1, type = "distribution")
predict(mm02, newdata = nd[1:5,], q = -2:2, 
        margins = 2, type = "distribution")

predict(mm1, newdata = nd[1:5,], q = -2:2, 
        margins = 1, type = "distribution")
predict(mm2, newdata = nd[1:5,], q = -2:2, 
        margins = 2, type = "distribution")

## expect correlations to be the same for the model with constant lambdas
c(coef(mm01, newdata = nd[1:5,], type = "Cor"))
c(coef(mm02, newdata = nd[1:5,], type = "Cor"))

## correlations for models with x-dependent lambda
c(coef(mm1, newdata = nd[1:5,], type = "Cor"))
c(coef(mm2, newdata = nd[1:5,], type = "Cor"))


##### mix of BoxCox and Colr margins: ##### 
d$Y1 <- (d$Y1 - min(d$Y1))/(max(d$Y1) - min(d$Y1))

b1 <- as.mlt(Colr(Y1 ~ X1 + X2 + X3, data = d, order = 1))
b2 <- as.mlt(Lm(Y2 ~ X1 + X2 + X3, data = d))

mm01 <- mmlt(b1, b2, formula = ~ 1, data = d)
mm02 <- mmlt(b2, b1, formula = ~ 1, data = d)

logLik(mm01)
logLik(mm02)

coef(b1)
coef(b2)
coef(mm01)
coef(mm02)
# remember that: lb <- (off <- 0.5) + X %*% (cf <- c(0, 2, 0))


## checking gradient
all.equal(c(numDeriv::grad(mm01$ll, coef(mm01))),c(mm01$sc(coef(mm01))), 
          check.attributes = FALSE, tol = 1e-4)
all.equal(c(numDeriv::grad(mm02$ll, coef(mm02))),c(mm02$sc(coef(mm02))), 
          check.attributes = FALSE, tol = 1e-4)

## covariate-dependent Lambda
mm1 <- mmlt(b1, b2, formula = ~ X1 + X2 + X3, data = d)
mm2 <- mmlt(b2, b1, formula = ~ X1 + X2 + X3, data = d)
logLik(mm1)
logLik(mm2)

coef(b1)
coef(b2)
coef(mm1)
coef(mm2)
# remember that: lb <- (off <- 0.5) + X %*% (cf <- c(0, 2, 0))

## checking gradient for diag = TRUE
all.equal(c(numDeriv::grad(mm1$ll, coef(mm1))),c(mm1$sc(coef(mm1))), 
          check.attributes = FALSE, tol = 1e-4)
all.equal(c(numDeriv::grad(mm2$ll, coef(mm2))),c(mm2$sc(coef(mm2))), 
          check.attributes = FALSE, tol = 1e-4)


N <- 1000
S <- diag(4)
S[lower.tri(S)] <- S[upper.tri(S)] <- .5

x <- matrix(runif(N*2), ncol = 2)

y <- x %*% matrix(c(1, -1, -.5, .5, -.2, .2, .3, -.3), nrow = 2) + rmvnorm(N, sigma = S)
d <- data.frame(y = y, x = x)

m1 <- Lm(y.1 ~ x.1 + x.2, data = d)
m2 <- Lm(y.2 ~ x.1 + x.2, data = d)
m3 <- Lm(y.3 ~ x.1 + x.2, data = d)
m4 <- Lm(y.4 ~ x.1 + x.2, data = d)

## simple formula
mc01 <- mmlt(m1, m2, m3, m4, formula = ~ 1, data = d)
mc02 <- mmlt(m2, m3, m1, m4, formula = ~ 1, data = d)

logLik(mc01)
logLik(mc02)

## complex formula
mc1 <- mmlt(m1, m2, m3, m4, formula = ~ x.1 + x.2, data = d)
mc2 <- mmlt(m2, m3, m1, m4, formula = ~ x.1 + x.2, data = d)

logLik(mc1)
logLik(mc2)

S <- diag(4)
x <- matrix(runif(N*2), ncol = 2)
S[lower.tri(S)] <- S[upper.tri(S)] <- .5

y <- x %*% matrix(c(1, -1, -.5, .5, -.2, .2, .3, -.3), nrow = 2) + rmvnorm(N, sigma = S)
d <- data.frame(y = y, x = x)

m1 <- Lm(y.1 ~ x.1 + x.2, data = d)
m2 <- Lm(y.2 ~ x.1 + x.2, data = d)
m3 <- Lm(y.3 ~ x.1 + x.2, data = d)
m4 <- Lm(y.4 ~ x.1 + x.2, data = d)
# m1$todistr$name <- m2$todistr$name <- m3$todistr$name <- m4$todistr$name <- "CF"

## simple formula
mc01 <- mmlt(m1, m2, m3, m4, formula = ~ 1, data = d, conditional = FALSE)

cf <- coef(mc01)
vr <- diag(vcov(mc01))
i <- grep("x", names(cf))

### same results
ret <- cbind(c(coef(m1), coef(m2), coef(m3), coef(m4)),
             cf[i],
             c(diag(vcov(m1)), diag(vcov(m2)), diag(vcov(m3)), diag(vcov(m4))),
             vr[i])
ret

#### check density

C <- as.array(coef(mc01, type = "Cor"))[,,1]
d1 <- sapply(1:N, function(i) dmvnorm(y[i,], mean = x[i,,drop = FALSE] %*% matrix(ret[,1], nrow = 2), sigma = C, log = TRUE))
d2 <- predict(mc01, newdata = d, type = "density", log = TRUE)
