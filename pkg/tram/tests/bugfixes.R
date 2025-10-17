
library("tram")
library("sandwich")
library("survival")
library("numDeriv")

### scores for fixed parameters
cars$int <- 1
a0 <- Lm(dist ~ speed, data = cars)
a1 <- Lm(dist ~ speed + int, data = cars, fixed = c("int" = 0))

s0 <- estfun(a0)
s1 <- estfun(a1)
s2 <- estfun(a1, parm = coef(as.mlt(a1), fixed = TRUE))

stopifnot(all.equal(s0[,"(Intercept)"], -s2[,"int"]))
stopifnot(all.equal(s1[,"(Intercept)"], -s2[,"int"]))
stopifnot(all.equal(s2[,"(Intercept)"], -s2[,"int"]))

### log_first for count data (by Sandra Siegfried)
### use: y + 1; log_first = TRUE, support = c(1, ...), bounds = c(1, ...)
set.seed(29)

Nsim <- 100
n <- 4000
b1 <- -4.5
b0 <- 5
theta <- 2

h <- qlogis(ppois(0:100, lambda = 5))

x <- runif(n, min = 0, max = 1) 
log.mu <- b0 + b1 * x

h.m <- matrix(h, nrow = length(h), ncol = length(x)) 
p <- (plogis(t(h.m) - b1 * x) - runif(n))^2
y <- max.col(-p) - 1

d <- data.frame(x = x, y = as.integer(y))
d$y.p1 <- d$y + 1L

m1 <- Colr(y ~ x, data = d,
                  support = c(0L, as.numeric(max(d$y))),
                  bounds = c(0L, Inf),
                  order = 10)

try(m2a <- Colr(y ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y))),
                  bounds = c(1L, Inf),
                  log_first = TRUE,
                  order = 10))

try(m2b <- Colr(y ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y))),
                  bounds = c(0L, Inf),
                  log_first = TRUE,
                  order = 10))

m3 <- Colr(y.p1 ~ x, data = d,
                  support = c(1L, as.numeric(max(d$y.p1))),
                  bounds = c(1L, Inf),
                  log_first = TRUE,
                  order = 10)

l1 <- m1$logliki(coef(as.mlt(m1)), rep(1, nrow(d)))
l3 <- m3$logliki(coef(as.mlt(m3)), rep(1, nrow(d)))
stopifnot(cor(l1, l3) > .9)

### 0.3-0
data("GBSG2", package = "TH.data")
### this model included an additional intercept
m1 <- Coxph(Surv(time, cens) | menostat:tgrade ~ horTh, data = GBSG2)
m2 <- Coxph(Surv(time, cens) | 0 + menostat:tgrade ~ horTh, data = GBSG2)
stopifnot(max(abs(coef(as.mlt(m1)) - coef(as.mlt(m2)))) < 
          sqrt(.Machine$double.eps))
### interaction term with stratum
m3 <- Coxph(Surv(time, cens) | horTh ~ menostat + menostat:horTh, 
            data = GBSG2)
ci <- confint(m3)["menostatPre:horThyes",]
stopifnot(all(!is.finite(ci)))

### problems with responses of class "R", spotted by Balint Tamasi
data("wine", package = "ordinal")
erating <- wine$rating
lrating <- erating
rrating <- erating
l9 <- lrating[wine$judge == 9]
l9[l9 > 1] <- levels(l9)[unclass(l9[l9 > 1]) - 1]
r9 <- rrating[wine$judge == 9]
r9[r9 < 5] <- levels(r9)[unclass(r9[r9 < 5]) + 1]
lrating[wine$judge != 9] <- rrating[wine$judge != 9] <- NA
erating[wine$judge == 9] <- NA
lrating[wine$judge == 9] <- l9
rrating[wine$judge == 9] <- r9
which(wine$judge == 9)
wine$crating <- R(erating, cleft = lrating, cright = rrating)
### gave an error
m <- Polr(crating ~ temp, data = wine, method = "probit")

### another one
data("GBSG2", package = "TH.data")
tmp <- GBSG2
tmp$y <- R(with(GBSG2, Surv(time, cens)))
m1 <- Coxph(Surv(time, cens) ~ horTh, data = tmp)
m2 <- Coxph(y ~ horTh, data = tmp)
stopifnot(all.equal(coef(as.mlt(m1)), coef(as.mlt(m2)), tol = 1e-4,
                    check.attributes = FALSE))

### check probabilistic index <-> log-odds ratio conversion
stopifnot(max(abs(PI(prob = PI(-15:15)) - (-15:15))) < 1e5)

### check updating with permutations
data("BostonHousing2", package = "mlbench")
m <- Colr(Surv(cmedv, cmedv < 50) ~ chas + crim, data = BostonHousing2)
mm <- as.mlt(m)
p <- sample(1:NROW(BostonHousing2))
## model with crim permuted
m1 <- update(mm, perm = "crim", permutation = p)
## permute data directly
tmp <- BostonHousing2
tmp[, "crim"] <- tmp[p, "crim"]
m2 <- Colr(Surv(cmedv, cmedv < 50) ~ chas + crim, data = tmp)
stopifnot(all.equal(coef(m1), coef(as.mlt(m2)), tol = 1e-3))
stopifnot(all.equal(logLik(m1), logLik(m2), tol = 1e-6))

### contraints, by Lucas Kook
data("GBSG2", package = "TH.data")
# gave an error
m <- Survreg(Surv(time, cens) ~ horTh + age, data = GBSG2, constraints = c("age >= 0"))

### mtram with interval censoring, spotted by Sandra Siegfried
dir <- system.file("rda", package = "TH.data")
load(file.path(dir, "Primary_endpoint_data.rda"))
trt <- "randarm5-FU + Oxaliplatin"
### convert "exact" event dates to interval-censoring (+/- one day)
tmp <- CAOsurv$iDFS
exact <- tmp[,3] == 1 
tmp[exact,2] <- tmp[exact,1] + 2
tmp[exact,1] <- pmax(tmp[exact,1] - 2, 0)
tmp[exact,3] <- 3
CAOsurv$iDFS2 <- tmp
CAO_SR <- Survreg(iDFS2 ~ randarm, data = CAOsurv, support = c(1, 1700), bounds = c(0, Inf))
CAO_SR_mtram <- mtram(CAO_SR, ~ (1 | Block), data = CAOsurv)
CAO_Cox <- Coxph(iDFS2 ~ randarm, data = CAOsurv, support = c(1, 1700), bounds = c(0, Inf), log_first = TRUE, order = 1)
CAO_Cox_mtram <- mtram(CAO_Cox, ~ (1 | Block), data = CAOsurv)
### wasn't equal
stopifnot(isTRUE(all.equal(logLik(CAO_SR_mtram), logLik(CAO_Cox_mtram), tol = 1e-5)))

## produces negative variances in shift-scale model
set.seed(100)
N <- 100
sc <- -5.5
sh <- 3
 
FZ <- pnorm
FZi <- qnorm
h2y <- function(y) log(-FZ(y, lower.tail = FALSE, log.p = TRUE))
 
x <- rep(x0 <- gl(2, 1), each = N)
xx <- (0:1)[x]
 
U <- runif(length(x))
scale <- sqrt(exp(sc * xx))
shift <- sh * xx
d <- data.frame(y = h2y( shift + FZi(U) / scale),
                x = x)
m <- BoxCox(y ~ x | x, data = d, scale_shift = TRUE)
coef(m)
try(diag(vcov(m)))

### constraints in remove intercept were incorrect
### spotted by Balint
stopifnot(isTRUE(all.equal(
  coef(as.mlt(Lm(dist ~ 1, data = cars, remove_intercept = FALSE)))[2:1] * c(-1, 1),
  coef(as.mlt(Lm(dist ~ 1, data = cars))), check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(
  coef(as.mlt(Survreg(dist ~ 1, data = cars, remove_intercept = FALSE)))[2:1] * c(-1, 1),
  coef(as.mlt(Survreg(dist ~ 1, data = cars))), check.attributes = FALSE)))

### Mmlt with shift-scale had incorrect gradients
m1 <- Lm(Sepal.Width ~ Petal.Length | Petal.Width, data = iris)
m2 <- Lm(Sepal.Length ~ Petal.Length | Petal.Width, data = iris)
m0 <- Mmlt(m1, m2, data = iris, formula = ~ 1, domargins = FALSE)
cf <- coef(m0)
m <- Mmlt(m1, m2, data = iris, formula = ~ 1, dofit = FALSE)
stopifnot(isTRUE(all.equal(c(-colSums(estfun(m, parm = cf))), c(grad(m$ll, cf)), check.attributes = FALSE)))

m1 <- Lm(Sepal.Width ~ Petal.Length | Petal.Width, data = iris, scale_shift = TRUE)
m2 <- Lm(Sepal.Length ~ Petal.Length | Petal.Width, data = iris, scale_shift = TRUE)
m0 <- Mmlt(m1, m2, data = iris, formula = ~ 1, domargins = FALSE)
cf <- coef(m0)
m <- Mmlt(m1, m2, data = iris, formula = ~ 1, dofit = FALSE)
stopifnot(isTRUE(all.equal(c(-colSums(estfun(m, parm = cf))), c(grad(m$ll, cf)), check.attributes = FALSE)))

### tram didn't allow binary factors
d <- data.frame(y = gl(2, 50), x = runif(100))
m0 <- tram(y ~ x, data = d, transformation = "discrete", distribution = "Logistic")
m1 <- Polr(y ~ x, data = d, method = "logistic")
m2 <- glm(y ~ x, data = d, family = binomial())
stopifnot(isTRUE(all.equal(c(logLik(m0)), c(logLik(m1)))))
stopifnot(isTRUE(all.equal(c(logLik(m0)), c(logLik(m2)))))

### requires mlt >= 1.5-3, which can handle missing values in the response
chk <- function(x, y, tol = 1e-3) stopifnot(all.equal(x, y, tol = tol, 
                                          check.attributes = FALSE))
N <- 50
d <- data.frame(y = rnorm(N), x = runif(N))
ic <- 1:10
d$y[ic] <- NA

m1 <- BoxCox(y ~ 1, data = d, na.action = na.pass)
m2 <- BoxCox(y ~ 1, data = d)
m3 <- BoxCox(y ~ 1, data = d[-ic,,drop = FALSE])
chk(logLik(m2), logLik(m1))
chk(logLik(m3), logLik(m1))
chk(nrow(m1$data), N)
chk(nrow(m2$data), sum(!is.na(d$y)))

m1 <- BoxCox(y ~ x, data = d, na.action = na.pass)
m2 <- BoxCox(y ~ x, data = d)
m3 <- BoxCox(y ~ x, data = d[-ic,,drop = FALSE])
chk(logLik(m2), logLik(m1))
chk(logLik(m3), logLik(m1))
chk(nrow(m1$data), N)
chk(nrow(m2$data), sum(!is.na(d$y)))
chk(coef(m2), coef(m1), tol = 1e-2)
chk(coef(m3), coef(m1), tol = 1e-2)

### allow to extract all relevant model matrix parts
## tram
m <- Coxph(Surv(time, cens) | horTh ~ age, data = GBSG2)
model.matrix(m, data = GBSG2[1:10,], what = "shifting")
model.matrix(m, data = GBSG2[1:10,], what = "interacting")
## stram
m <- Coxph(Surv(time, cens) | horTh ~ age | pnodes, data = GBSG2)
model.matrix(m, data = GBSG2[1:10,], what = "shifting")
model.matrix(m, data = GBSG2[1:10,], what = "scaling")
model.matrix(m, data = GBSG2[1:10,], what = "interacting")

### check numerically approximated Hessians when analytic Hessian is
### not available
N <- 100
x <- runif(N)
y <- rnorm(N)
d <- data.frame(x = x, y = y)

### no scaling, optimHess
m1 <- Colr(y ~ x | x, data = d, scaleparm = TRUE, 
  optim = mltoptim(hessian = FALSE))
### no scaling, Hessian from auglag
m2 <- Colr(y ~ x | x, data = d, scaleparm = TRUE, 
  optim = mltoptim(hessian = TRUE))
### scaling, optimHess
m3 <- Colr(y ~ x | x, data = d, scaleparm = FALSE, 
  optim = mltoptim(hessian = FALSE))
### scaling, Hessian from auglag
m4 <- Colr(y ~ x | x, data = d, scaleparm = FALSE, 
  optim = mltoptim(hessian = TRUE))

stopifnot(all.equal(logLik(m1), logLik(m2), logLik(m3), logLik(m4)))

cf0 <- coef(m1)
cf <- coef(as.mlt(m1))
i <- match(names(cf0), names(cf))

f <- function(x) logLik(m1, parm = x)
stopifnot(max(abs(solve(-hessian(f, cf))[i,i] - vcov(m1))) < 1e-5)

f <- function(x) logLik(m2, parm = x)
stopifnot(max(abs(solve(-hessian(f, cf))[i,i] - vcov(m2))) < 1e-5)

f <- function(x) logLik(m3, parm = x)
stopifnot(max(abs(solve(-hessian(f, cf))[i,i] - vcov(m3))) < 1e-3)

### optimHess and numDeriv::hessian disagree a bit
f <- function(x) logLik(m4, parm = x)
stopifnot(max(abs(solve(-hessian(f, cf))[i,i] - vcov(m4))) < 1e-1)
