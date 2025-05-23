
library("tram")
library("sandwich")
set.seed(29)

chk <- function(x, y, tol = 1e-3) 
    stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE, tol = tol)))

if (require("mlbench")) {

  data("BostonHousing2", package = "mlbench")

  m <- Lm(cmedv ~ nox, data = BostonHousing2)

  d <- predict(as.mlt(m), type = "density")
  ll <- sum(log(d))
  chk(c(logLik(m)), ll)
  d1 <- predict(as.mlt(m), newdata = BostonHousing2[1:6,], type = "density")
  d2 <- diag(predict(as.mlt(m), newdata = BostonHousing2[1:6,-6], q =
             BostonHousing2[1:6,6], type = "density"))
  chk(d1, d2)
  prb <- 1:9 / 10
  q <- predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], prob = prb,
               type = "quantile")
  p1 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,1], 
                      type = "distribution")[,1], 2)
  p2 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,2], 
                      type = "distribution")[,2], 2)
  chk(p1, prb)
  chk(p2, prb)

  m <- Lm(cmedv ~ nox | rm, data = BostonHousing2)

  d <- predict(as.mlt(m), type = "density")
  ll <- sum(log(d))
  chk(c(logLik(m)), ll)
  d1 <- predict(as.mlt(m), newdata = BostonHousing2[1:6,], type = "density")
  d2 <- diag(predict(as.mlt(m), newdata = BostonHousing2[1:6,-6], q =
             BostonHousing2[1:6,6], type = "density"))
  chk(d1, d2)
  prb <- 1:9 / 10
  q <- predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], prob = prb,
               type = "quantile")
  p1 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,1], 
                      type = "distribution")[,1], 2)
  p2 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,2], 
                      type = "distribution")[,2], 2)
  chk(p1, prb)
  chk(p2, prb)

}

if (suppressPackageStartupMessages(require("gamlss"))) {

  N <- 1000
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  y <- rnorm(N, mean = 2 * x1, sd = 2 * sqrt(exp(.5 * x2)))
  d <- data.frame(y = y, x1 = x1, x2 = x2)

  ### shift-scale transformation model
  tm <- Lm(y ~ x1 | x2, data = d, scale_shift = TRUE)

  ### same model, via gamlss
  gm <- gamlss(y ~ x1, sigma.formula = ~ x2, data = d, family = NO2, trace = FALSE)
  chk(logLik(gm), logLik(tm))
  chk(-coef(tm)["scl_x2"], coefAll(gm)[[2]]["x2"])

  X1 <- model.matrix(~ x1, data = d)
  X2 <- model.matrix(~ x2, data = d)

  #### logLik gamlss
  glp1 <- X1 %*% coefAll(gm)[[1]]
  glp2 <- X2 %*% coefAll(gm)[[2]]
  chk(c(logLik(gm)), 
      sum(dnorm(y, mean = glp1, sd = sqrt(exp(glp2)), log = TRUE)))

  ### shift-scale transformation model
  Y <- model.matrix(~ y, data = d)
  h0 <- Y %*% coef(as.mlt(tm))[1:2]
  cf <- coef(tm)
  tlp1 <- X1[,-1,drop = FALSE] %*% cf[-grep("scl", names(cf))]
  tlp2 <- X2[,-1,drop = FALSE] %*% cf[grep("scl", names(cf))]
  st <- exp(tlp2)
  xi <- coef(as.mlt(tm))["y"]

  ### transformation function
  chk(c(sqrt(st) * (h0 - tlp1)), tr <- predict(tm, type = "trafo"))

  ### logLik
  chk(c(logLik(tm)), 
      sum(log(1 / sqrt(2 * pi)) + .5 * tlp2 + log(xi) - 1 / 2 * tr^2))

  ### score function
  # score alpha
  sa <- sum(-sqrt(st) * tr)
  # score xi
  sx <- sum(1 / xi - sqrt(st) * tr * y)
  # score beta
  sb <- sum(sqrt(st) * tr * X1[,-1,drop = FALSE])
  # score gamma
  sg <- sum(1 / 2 * (1 - tr^2) * X2[,-1, drop = FALSE])

  chk(c(sa, sx, sb, sg), -Gradient(tm))

  ### Hessian
  # alpha, alpha
  Haa <- sum(-st)
  # xi, xi
  Hxx <- sum(- xi^(-2) - st * y^2)
  # beta, beta
  Htm <- crossprod(-st * X1[, -1, drop = FALSE], X1[, -1, drop = FALSE])
  # gamma, gamma
  Hgg <- crossprod(-1 / 2 * tr^2 * X2[, -1, drop = FALSE], X2[, -1, drop = FALSE])
  # alpha, xi
  Hax <- sum(-st * y)
  # alpha, beta
  Hab <- colSums(st * X1[, -1, drop = FALSE])
  # alpha, gamma
  Hag <- colSums(- sqrt(st) * tr * X2[, -1, drop = FALSE])
  # beta, xi
  Hbx <- colSums(st * y * X1[, -1, drop = FALSE])
  # beta, gamma
  Hbg <- crossprod(sqrt(st) * tr * X1[, -1, drop = FALSE], X2[, -1, drop = FALSE])
  # xi, gamma
  Hxg <- colSums(- sqrt(st) * tr * y * X2[, -1, drop = FALSE])
  H <- diag(c(Haa, Hxx, Htm, Hgg))
  H[1, 2] <- H[2, 1] <- Hax
  H[1, 3] <- H[3, 1] <- Hab
  H[1, 4] <- H[4, 1] <- Hag
  H[2, 3] <- H[3, 2] <- Hbx
  H[2, 4] <- H[4, 2] <- Hxg
  H[3, 4] <- H[4, 3] <- Hbg

  chk(-H, Hessian(as.mlt(tm)))

  ### scale terms are defined analoguously in Lm and gamlss(NO2)
  ### compare z statistics
  chk(-coef(tm)["scl_x2"] / sqrt(vcov(tm)["scl_x2", "scl_x2"]),
      coefAll(gm)[[2]]["x2"] / sqrt(vcov(gm)["x2", "x2"]))
}

### predict

N <- 10000
x1 <- runif(N)
x2 <- runif(N)
s <- gl(k <- 11, N/k)
y <- rnorm(N, mean = 2 * x1, sd = sqrt(exp(.5 * x2)))
d <- data.frame(y = y, x1 = x1, x2 = x2, s = s)

fm <- formula(y ~ s + x1)
m <- Coxph(fm, data = d)
try(a <- predict(m, type = "distribution")) ### works

sfm <- as.formula(y ~ s | x1)
sm <- Coxph(sfm, data = d)
try(b <- predict(sm, type = "distribution")) ### didn't work

### check discrete models
data("wine", package = "ordinal")
m <- Polr(rating ~ temp | contact, data = wine, 
          optim = mltoptim(spg = list(maxit = 10000, checkGrad = TRUE)))

### inference methods & residuals
N <- 1000
x1 <- runif(N)
x2 <- runif(N)
z <- rnorm(N)
y <- (z + .2 * x1) * sqrt(exp(.5 * x2))
d <- data.frame(y = y, x1 = x1, x2 = x2, one = 1)

### Wald, Score and Likelihood intervals
m2 <- Lm(y ~ x1 | x2, data = d)
ciW <- confint(m2)
ciL <- confint(profile(m2))
ciS <- rbind(score_test(m2, "x1")$conf,
             score_test(m2, "scl_x2")$conf)

chk(ciW, ciL, tol = 1e-3)
chk(ciW, ciS, tol = 1e-3)

### not quite the same!
round(ciW[1,], 3)
round(perm_test(m2, "x1")$conf, 3)

### check residuals
m0 <- Lm(y ~ 1, data = d)
rs <- c(estfun(as.mlt(m0)) %*% coef(as.mlt(m0))) * .5
rs2 <- resid(m0, what = "scaling")
chk(rs, rs2)

### residuals, nonparametrically
N <- 50
x1 <- runif(N)
x2 <- runif(N)
z <- rnorm(N)
y <- (z + .2 * x1) * sqrt(exp(.5 * x2))
d <- data.frame(y = y, x1 = x1, x2 = x2, one = 1)

d$yR <- R(y, as.R.ordered = TRUE)
m0 <- Polr(yR ~ 1, data = d)
rs <- c(estfun(as.mlt(m0)) %*% coef(as.mlt(m0))) * .5
rs2 <- resid(m0, what = "scaling")
chk(rs, rs2)

### linear predictor
aq <- airquality[complete.cases(airquality),]
aq <- as.data.frame(lapply(aq, as.double))
m <- Lm(Ozone ~ Solar.R + Wind + Temp + Month + Day | Solar.R + Wind + Temp + Month + Day, 
        data = aq)

lp1 <- predict(m, type = "lp", what = "shifting")
lp2 <- predict(m, type = "lp", what = "scaling")

X <- -model.matrix(m)
stopifnot(identical(lp1, X %*% coef(m, with_baseline = FALSE)[m$shiftcoef]))
X <- model.matrix(m, what = "scaling")
stopifnot(identical(lp2, X %*% coef(m, with_baseline = FALSE)[m$scalecoef]))
