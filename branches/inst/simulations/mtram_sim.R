### mtram simulations

## packages
library("geepack")
library("lme4")
library("mvtnorm")
library("tram")
library("tramME")


## some functions
PFUN <- function(x, ...) pchisq(x, df = 9, ...)
QFUN <- function(p, ...) qchisq(p, df = 9, ...)
DFUN <- function(x, ...) dchisq(x, df = 9, ...)
h1 <- function(x) QFUN(plogis(x, log.p = TRUE), log.p = TRUE) 
#                 ^^^^^^^^^ <- h^{-1}

## data generating process
## needs a bunch of variables to be defined globally!
dgp <- function(dichotomize = FALSE,
                tau, S, Sigma, D, fct) {
  
  Z <- rmvnorm(N, sigma = Sigma)
  x <- matrix(runif(N * Ni * p), ncol = p)
  h1 <- function(x) QFUN(plogis(x, log.p = TRUE), log.p = TRUE) 
  #                 ^^^^^^^^^ <- h^{-1}
  y <- h1(c(D %*% qlogis(pnorm(solve(D) %*% t(Z)))) / fct + x %*% beta)
  if(dichotomize) {
    ret <- data.frame(y = factor(as.integer(y > median(y))), x, cls = cls)
  } else {
    ret <- data.frame(y = y, x, cls = cls)
  }
  return(ret)
}

## model fitting wrapped in a function to catch errors
fits1 <- function(i) {
  
  d <- dgp(dichotomize = TRUE, tau, S, Sigma, D, fct)
  
  cm <- ctm(as.basis(~ y, data = d), shifting = ~ X1 + X2 + X3,
            data = d, todistr = "Logistic", negative = TRUE)
  m <- mlt(cm, data = d, fixed = c("y1" = 0))
  mt <- mtram(m, ~ (1 | cls), data = d, Hessian = TRUE)
  
  cf <- coef(mt)
  
  gm <- geeglm(I((0:1)[y]) ~ X1 + X2 + X3, id = cls, data = d, family = binomial,
               corstr = "exchangeable")
  
  # cfm[i,] <- cf[c("X1", "X2", "X3")] / sqrt(1 + cf["gamma1"]^2)
  # cfgee[i,] <- coef(gm)[-1L]

  ### compute confidence intervals for marginal odds ratios from mtram
  Z <- rmvnorm(10000, mean = coef(mt)[-1L], sigma = solve(mt$Hessian)[-1,-1])
  ci <- apply(Z[,-4] / sqrt(1 + Z[,4]^2), 2, quantile, prob = c(.025, .975))
  
  return(list(cfm = cf[c("X1", "X2", "X3")] / sqrt(1 + cf["gamma1"]^2),
              cfgee = coef(gm)[-1L],
              lwrm = ci[1,], uprm = ci[2,],
              lwrgee = coef(gm)[-1L] - qnorm(.975) * summary(gm)[["coefficients"]][-1L, "Std.err"],
              uprgee = coef(gm)[-1L] + qnorm(.975) * summary(gm)[["coefficients"]][-1L, "Std.err"]))
  
}


fits2 <- function(i) {
  
  d <- dgp(dichotomize = FALSE, tau, S, Sigma, D, fct)
  
  m <- Colr(y ~ X1 + X2 + X3, data = d)
  mt <- mtram(m, ~ (1 | cls), data = d, Hessian = TRUE)
  
  cf <- coef(mt)
  
  Z <- rmvnorm(10000, mean = coef(mt)[c("X1", "X2", "X3", "gamma1")], 
               sigma = solve(mt$Hessian)[-(1:7),-(1:7)])
  ci <- apply(Z[,-4] / sqrt(1 + Z[,4]^2), 2, quantile, prob = c(.025, .975))
  
  return(list(cfm = cf[c("X1", "X2", "X3")] / sqrt(1 + cf["gamma1"]^2),
              lwrm = ci[1,], uprm = ci[2,]))
}


##### comparison with GEE #####
### marginal binary logit models

set.seed(150594)
Nsim <- 10000

## simulation parameters
N <- 100 # 100 clusters of
Ni <- 5 # size 5
cls <- gl(N, Ni)
p <- 3
beta <- c(0, 1, 2)
taus <- c(0, 0.5, 1, 1.5, 2, 3)

results <- vector(mode = "list", length = length(taus))



## on the server, the seed is reset for every new value of tau
for(j in 1:length(taus)) {
  
  set.seed(150594)
  tau <- taus[j]
  print(tau)
  
  cfm <- cfgee <- lwrm <- uprm <- lwrgee <- uprgee <- matrix(NA, ncol = p, nrow = Nsim)
  
  Ui <- matrix(1, ncol = 1, nrow = Ni)
  S <- tau^2
  Sigma <- S * tcrossprod(Ui) + diag(Ni)
  D <- diag(sqrt(diag(Sigma)))
  fct <- sqrt(1 + tau^2)
  
  for (i in 1:Nsim) {
    
    rets <- try(fits1(i), TRUE)
    while(inherits(rets, 'try-error')){
      rets <- try(fits1(i), TRUE)
    }
    print(i)
    
    cfm[i,] <- rets$cfm
    cfgee[i,] <- rets$cfgee
   
    lwrm[i,] <- rets$lwrm
    uprm[i,] <- rets$uprm
    
    lwrgee[i,] <- rets$lwrgee
    uprgee[i,] <- rets$uprgee
  }
  
  results[[j]]$cfm <- cfm
  results[[j]]$lwrm <- lwrm
  results[[j]]$uprm <- uprm
  results[[j]]$biasm <- rowMeans((t(cfm) - beta)^2)
  results[[j]]$CIlengthm <- colMeans(uprm - lwrm)
  results[[j]]$coveragem <- rowMeans(t(lwrm) < beta & beta < t(uprm))
  
  results[[j]]$cfgee <- cfgee
  results[[j]]$lwrgee <- lwrgee
  results[[j]]$uprgee <- uprgee
  results[[j]]$biasgee <- rowMeans((t(cfgee) - beta)^2)
  results[[j]]$CIlengthgee <- colMeans(uprgee - lwrgee)
  results[[j]]$coveragegee <- rowMeans(t(lwrgee) < beta & beta < t(uprgee))
}

dump(c("results"), file = "mtram_gee_dichotomized.R")

## extracting relevant info
sapply(results, function (x) c(sd(x$cfm[, 1]), sd(x$cfm[, 2]), sd(x$cfm[, 3])))
sapply(results, function (x) c(sd(x$cfgee[, 1]), sd(x$cfgee[, 2]), sd(x$cfgee[, 3])))

sapply(results, function (x) x$biasm)
sapply(results, function (x) x$biasgee)

sapply(results, function (x) x$CIlengthm)
sapply(results, function (x) x$CIlengthgee)

sapply(results, function (x) x$coveragem)
sapply(results, function (x) x$coveragegee)



##### mtram with non-dichotomized data #####

set.seed(150594)
Nsim <- 10#000

## simulation parameters
N <- 100 # 100 clusters of
Ni <- 5 # size 5
cls <- gl(N, Ni)
p <- 3
beta <- c(0, 1, 2)
taus <- c(0, 0.5, 1, 1.5, 2, 3)

results <- vector(mode = "list", length = length(taus))


for(j in 1:length(taus)) {
  
  set.seed(150594)
  tau <- taus[j]
  print(tau)
  
  cfm <- lwrm <- uprm <- matrix(NA, ncol = p, nrow = Nsim)
  
  Ui <- matrix(1, ncol = 1, nrow = Ni)
  S <- tau^2
  Sigma <- S * tcrossprod(Ui) + diag(Ni)
  D <- diag(sqrt(diag(Sigma)))
  fct <- sqrt(1 + tau^2)
  
  for (i in 1:Nsim) {
    
    rets <- try(fits2(i), TRUE)
    while(inherits(rets, 'try-error')){
      rets <- try(fits2(i), TRUE)
    }
    print(i)
    
    cfm[i,] <- rets$cfm
    
    lwrm[i,] <- rets$lwrm
    uprm[i,] <- rets$uprm
  }
  
  results[[j]]$cfm <- cfm
  results[[j]]$lwrm <- lwrm
  results[[j]]$uprm <- uprm
  results[[j]]$bias <- rowMeans((-t(cfm) - beta)^2)
  results[[j]]$CIlength <- colMeans(uprm - lwrm)
  results[[j]]$coverage <- rowMeans(t(lwrm) < -beta & -beta < t(uprm))
  
}

dump(c("results"), file = "mtram_continuous.R")


## extracting relevant info
sapply(results, function (x) c(sd(x$cfm[, 1]), sd(x$cfm[, 2]), sd(x$cfm[, 3])))

sapply(results, function (x) x$bias)

sapply(results, function (x) x$CIlength)

sapply(results, function (x) x$coverage)


##### comparison with tramME #####
### compare mtram and tramME in a marginally interpretable for h = qlogis chisq
### model, expect mtram to outperform tramME

## --- Numerical integration
## A function to evaluate the joint cdf of the response and the random effects:
## Takes a vector of random effect and covariates values, evaluates the conditional
## distribution at these values and multiplies it with the pdf of the random effects
jointCDF <- function(re, nd, mod, mod_fe) {
  nd <- nd[rep(1, length(re)), ]
  nd$cls <- seq(nrow(nd)) ## to take vector-valued REs
  pr <- predict(mod, newdata = nd, ranef = re, type = "distribution") *
    dnorm(re, 0, sd = sqrt(varcov(mod)[[1]][1, 1]))
  c(pr)
}
## Marginalize the joint cdf by integrating out the random effects
## using adaptive quadratures
marginalCDF <- function(nd, mod, mod_fe) {
  
  if(sqrt(varcov(mod)[[1]][1, 1]) > .05) {
    nd$cdf <- integrate(jointCDF, lower = -Inf, upper = Inf, nd = nd, 
                      mod = mod, mod_fe = mod_fe)$value
  } else { ## predict from unconditional model to get sensible predictions
    mod_fe$coef <- mod$param$beta[1:length(mod_fe$coef)]
    nd$cdf <- predict(mod_fe, newdata = nd, type = "distribution")
  }
  nd
}
MC <- 4L


set.seed(150594)
Nsim <- 100

## simulation parameters
N <- 100 # 100 clusters of
Ni <- 5 # size 5
cls <- gl(N, Ni)
p <- 3
beta <- c(0, 1, 2)
taus <- c(0, 0.5, 1, 1.5, 2, 3)

mm_fe <- mBC_fe <- list()

f <- function(tau) {
  
  print(tau)
  
  Ui <- matrix(1, ncol = 1, nrow = Ni)
  S <- tau^2
  Sigma <- S * tcrossprod(Ui) + diag(Ni)
  D <- diag(sqrt(diag(Sigma)))
  fct <- sqrt(1 + tau^2)
  
  d <- dgp(dichotomize = FALSE, tau, S, Sigma, D, fct)
  m <- Colr(y ~ X1 + X2 + X3, data = d)
  mtT <- mtram(m, ~ (1 | cls), data = d, Hessian = FALSE)
  mm <- ColrME(y ~ X1 + X2 + X3 + (1 | cls), data = d)
  mm_fe <- Colr(y ~ X1 + X2 + X3, data = d)
  mBC <- BoxCoxME(y ~ X1 + X2 + X3 + (1 | cls), data = d, order = 1)
  mBC_fe <- BoxCox(y ~ X1 + X2 + X3, data = d, order = 1)
  
  ### quantiles true marginal distribution
  x <- rep(.5, 3)
  p <- 1:99 / 100
  iq <- h1(qlogis(p) + c(x %*% beta))
  y <- d$y
  h <- function(y) qlogis(PFUN(y))
  hp <- function(y) dlogis(PFUN(y)) * DFUN(y)
  d <- dlogis(h(iq) - c(x %*% beta)) * hp(iq)
  
  ## Set up the grid on which we evaluate the marginal distribution
  nd <- expand.grid(y = iq,
                    X1 = 0.5, X2 = 0.5, X3 = 0.5)
  nd$d <- d
  nd$p <- p
  nd1 <- nd
  ## Calls marginalCDF on each row of nd
  ## (done in parallel to speed up computations)
  mp_mm <- parallel::mclapply(split(nd1, seq(nrow(nd1))),
                              marginalCDF, mod = mm, mod_fe = mm_fe, 
                              mc.cores = MC)
  mp_mm <- do.call("rbind", mp_mm)
  
  mp_mBC <- parallel::mclapply(split(nd1, seq(nrow(nd1))),
                               marginalCDF, mod = mBC, mod_fe = mBC_fe,
                               mc.cores = MC)
  mp_mBC <- do.call("rbind", mp_mBC)
  
  nd <- merge(nd, mp_mm[, c("y", "cdf")])
  nd$cdfBC <- mp_mBC$cdf
  
  ## mtram
  m1 <- m
  m1$coef <- coef(mtT)[-11]
  m1$coef[1:10] <- m1$coef[1:10]/sqrt(1 + coef(mtT)[11]^2)
  nd$cdfmtram <- predict(m1, newdata = nd, type = "distribution")
  
  return(nd)
  
}


results <- vector(mode = "list", length = Nsim)
for(j in 1:Nsim) {
  print(j)
  results[[j]] <- lapply(taus, f)
}

# dump(c("results"), file = "res_mtram_tramME.R")
save(results, file = "res_mtram_tramME.RData")
# load("res_mtram_tramME.RData")

### plots with lattice
library("lattice")
library("latticeExtra")
par(las = 1)

## difference in cdfs
res_l <- lapply(results, function(x) {
  
  nd <- rbind(x[[1]][, 1:6], x[[1]][, 1:6], x[[1]][, 1:6])
  nd$cdf <- c(x[[1]]$cdf, x[[1]]$cdfBC, x[[1]]$cdfmtram)
  nd$model <- factor(c(rep("ColrME", 99), rep("lmer", 99), rep("mtram", 99)))
  nd$tau <- factor(c(rep(taus[1], 99*3)))

  for(i in 2:length(taus)) {
    nd1 <- rbind(x[[i]][, 1:6], x[[i]][, 1:6], x[[i]][, 1:6])
    nd1$cdf <- c(x[[i]]$cdf, x[[i]]$cdfBC, x[[i]]$cdfmtram)
    nd1$model <- factor(c(rep("ColrME", 99), rep("lmer", 99), rep("mtram", 99)))
    nd1$tau <- factor(c(rep(taus[i], 99*3)))

    nd <- rbind(nd, nd1)
  }
  return(nd)
})

nd_all <- res_l[[1]]
for (i in 2:length(res_l)) {
  nd_all <- rbind(nd_all, res_l[[i]])
}

panel_f1 <- function(x, y, repl = 100) {
  yseq <- unique(res_l[[1]]$y)
  # panel.grid(h = -1, v = -1)
  for (i in 1:(repl+1)) {
    idx <- ((1+length(yseq)*(i-1)):(length(yseq)*i))
    panel.lines(x[idx], y[idx], col = "black", alpha = .1)
  }
  panel.abline(h = 0, col = "red", lty = 2)
}

# pdf(file = "sim-cdfdiff.pdf", height = 4.5)
xyplot(p - cdf ~ y | model + tau, data = nd_all,
       type = "l", 
       # between = list(x = 1), 
       layout = c(6,3),
       col = rgb(.1, .1, .1, 0),
       panel = panel_f1,
       xlab = "y", ylab = c(expression(paste("F(y) - ", hat(F), "(y)", sep = ""))),
       ylim = c(-.2, .2),
       strip = strip.custom(bg = "transparent")
       # strip = strip.custom(which.given = 2,
       #                      # var.name = "tau",
       #                      bg = "transparent",
       #                      factor.levels = c(expression(paste(gamma[1], " = 0", sep = "")),
       #                                        expression(paste(gamma[1], " = 0.5", sep = "")),
       #                                        expression(paste(gamma[1], " = 1", sep = "")),
       #                                        expression(paste(gamma[1], " = 1.5", sep = "")),
       #                                        expression(paste(gamma[1], " = 2", sep = "")),
       #                                        expression(paste(gamma[1], " = 3", sep = "")),
       #                                        "lmer", "mtram", "ColrME")),
       # strip = strip.custom(which.given = 1, bg = "transparent",
       #                      factor.levels = rep(levels(nd_all$model), 6)),
)
# dev.off()

## integrated MSE for mtram and ColrME
iMSE_s <- vector(mode = "list", length = length(taus))

for (i in 1:length(taus)) {
  
  iMSE <- lapply(results, function(y) {
    x <- y[[i]]
    c(integrate(splinefun(x$y, (x$cdf - x$p)^2 * x$d), lower = 0, upper = max(x$y))$value,
      integrate(splinefun(x$y, (x$cdfmtram - x$p)^2 * x$d), lower = 0, upper = max(x$y))$value)
  })
  iMSE_s[[i]] <- (matrix(unlist(iMSE), ncol = 2, byrow = TRUE))
}

nd_iMSE <- data.frame(ColrME = c(sapply(iMSE_s, function(x) {return(x[, 1])})),
                      mtram = c(sapply(iMSE_s, function(x) {return(x[, 2])})),
                      tau = factor(c(rep(taus, each = 100))))

panel_f2 <- function(x, y, ...) {
  panel.abline(a = 0, b = 1, col = "red", lty = 2)
  panel.xyplot(x, y, ...)
}

rng <- round(c(max(nd_iMSE$ColrME), max(nd_iMSE$mtram)), 4)

# pdf(file = "sim-imse.pdf", height = 4.5)
xyplot(ColrME ~ mtram | tau, data = nd_iMSE, 
       xlim = c(0, max(rng)), ylim = c(0, max(rng)),
       pch = 20, col = rgb(.1, .1, .1, .2),
       main = "Integrated MSE",
       xlab = "mtram", ylab = "ColrME",
       layout = c(3, 2),
       panel = panel_f2,
       strip = strip.custom(bg = "transparent",
                            factor.levels = c(expression(paste(gamma[1], " = 0", sep = "")),
                                              expression(paste(gamma[1], " = 0.5", sep = "")),
                                              expression(paste(gamma[1], " = 1", sep = "")),
                                              expression(paste(gamma[1], " = 1.5", sep = "")),
                                              expression(paste(gamma[1], " = 2", sep = "")),
                                              expression(paste(gamma[1], " = 3", sep = "")))),
)
# dev.off()
