### Code from
### "Smooth Transformation Models for Survival Analysis: A Tutorial Using R"
###   by Sandra Siegfried, Balint Tamasi & Torsten Hothorn

# required packages
pkgs <- c("mlt", "tram",  "trtf", "SparseGrid", "ATR", "tramME", "multcomp",
  "coin", "TH.data", "survival", "colorspace", "xtable")

ix <- which(!sapply(pkgs, require, char = TRUE))
if (length(ix) > 0) {install.packages(pkgs[ix], repos = "https://stat.ethz.ch/CRAN/")
 sapply(pkgs[ix], require, char = TRUE)}

set.seed(290875)

## plotting
xlab <- "Time (in days)"
lxlab <- paste0(xlab, " on log-scale")
ylabS <- "Probability of survival"
ylablHaz <- "Log-cumulative hazard"
ylabcumHR <- expression(Lambda[1](t)/Lambda[0](t))
ylimS <- c(0, 1)
ylimHR <- c(0, 1.6)
q <- 0:2204
xlim <- c(0, max(q))
lwd <- 1.3

## color
acol <- sequential_hcl(6, "BluYl")[1:5]
col <- acol[c(2, (length(acol)) - 1)]
lcol <- lighten(col, amount = .4) ## lighten color for overlaid lines 

## aux
perm_test_biv.stram <-  function(object, seed = 1) {
  stopifnot(inherits(object, "stram"))
  fixed <- c(trt = 0, scl = 0)
  lhs <- object$call[[2]][[3]]
  if (!(length(lhs) == 3 & lhs[[2]] == lhs[[3]]))
    stop("Bivariate score perm test not applicable")
  names(fixed) <- names(coef(object))
  m0 <- update(object, fixed = fixed) ## uncond. model
  r <- resid(m0, what = "shifting")
  rs <- resid(m0, what = "scaling")
  set.seed(seed)
  
  formula <- as.formula(paste("r + rs ~", lhs[[2]]))
  pvalue(independence_test(formula, data = m0$data))
}

## formatting
big.mark <- "'"
frmt0 <- round
frmt <- function(digits, x, math = FALSE) {
  if (!is.numeric(x)) return(x)
    ret <- formatC(round(x, digits), digits = digits, format = "f", big.mark = big.mark) 
    if (math) ret <- paste("$", ret, "$")
    if (is.matrix(x)) {
        ret <- matrix(ret, nrow = nrow(x))
        dimnames(ret) <- dimnames(x)
    }
    ret
}

frmt1 <- function(x, math = FALSE) frmt(1, x = x, math = math)
frmt2 <- function(x, math = FALSE) frmt(2, x = x, math = math)
frmt3 <- function(x, math = FALSE) frmt(3, x = x, math = math)

## logLik
frmtll <- function(x, math = FALSE, mark = FALSE) {
  if (!inherits(x, "logLik") && !is.numeric(x) && all(!is.na(x))) x <- logLik(x)
    if (is.na(x)) return("")
  ret <- frmt2(abs(x), math = FALSE)
  if (x < 0) ret <- paste0(ifelse(math, "$-$", "-"), ret)
  if (mark) ret <- paste0("{\\color{darkgray}", ret, "}")
  ret
}

## data
load(system.file("rda", "Primary_endpoint_data.rda", package = "TH.data"))

## randomization arm
levs <- levels(CAOsurv$randarm)
trt <- with(CAOsurv, paste0("randarm", levs[2], collapse = ""))
nd1 <- data.frame(randarm = factor(levs, levels = levs))

## strata
CAOsurv$strat <- with(CAOsurv, interaction(strat_t, strat_n))
slevs <- levels(CAOsurv$strat)
nd2 <- data.frame(randarm = nd1$randarm[1], strat = factor(slevs, levels = slevs))

## for pretty legends
lslevs <- gsub("\\.", " : ", slevs)
lslevs <- gsub("cT4", "cT4    ", lslevs)

## id
CAOsurv$id <- factor(1:nrow(CAOsurv))

## ----pars, include = FALSE----------------------------------------------------
par_main <- expression(par(mgp = c(2.5, 1, 0), mar = c(4, 4, 1.5, 4), las = 1))
par_surv <- expression(par(mgp = c(2.5, 1, 0), mar = c(6, 6, .5, 4), las = 1))



## ----packages, echo = FALSE---------------------------------------------------
library("tram")

## ----risk-tab-----------------------------------------------------------------
risktab <- function(ti, st) { ## time-index and survival time
  nrisk <- NULL
  for (t in ti) nrisk <- c(nrisk, sum(st >= t))
  return(nrisk)
}

plot.risktab <- function(tvar, ti = seq(min(q), max(q), by = 500),
  cex = .8, at = -450) {
mtext(levs[1], 1, line = 4, at = at, cex = cex)
mtext(risktab(ti, CAOsurv[CAOsurv$randarm == levs[1], tvar]),
  side = 1, line = 4, at = ti, cex = cex)
mtext(levs[2], 1, line = 5, at = at, cex = cex)
mtext(risktab(ti, CAOsurv[CAOsurv$randarm == levs[2], tvar]),
  side = 1, line = 5, at = ti, cex = cex)
}

## ----surv-OS------------------------------------------------------------------
surv_OS <- survfit(OS ~ randarm, data = CAOsurv) ## KM


## ----surv-iDFS----------------------------------------------------------------
surv_iDFS <- survfit(iDFS ~ randarm, data = CAOsurv) ## Turnbull






## ----CAO-table, results = 'hide'----------------------------------------------
tab <- xtabs( ~ strat + randarm, data = CAOsurv)
tab <- rbind(tab, "Total" = colSums(tab))



## ----surv-iDFS-plot-----------------------------------------------------------
eval(par_surv)
plot(surv_iDFS, ylim = ylimS, xlim = xlim,
  col = lcol, lwd = lwd, xlab = xlab, ylab = ylabS)
legend("bottomright", legend = levs, col = col, bty = "n", lty = 1, lwd = 1, cex = .8)
plot.risktab(tvar = "iDFStime")


## ----surv-OS-plot-------------------------------------------------------------
plot(surv_OS, ylim = ylimS, xlim = xlim,
  col = lcol, lwd = lwd, xlab = xlab, ylab = ylabS)
legend("bottomright", legend = levs, col = col, bty = "n", lty = 1, lwd = 1, cex = .8)
plot.risktab(tvar = "OStime")



## ----WEI-model-fit, echo = FALSE, cache = TRUE--------------------------------
mw <- 
Survreg(iDFS ~ randarm, data = CAOsurv, dist = "weibull")


## ----WEI-summary, cache = TRUE, results = "hide", fig.show = 'hide'-----------
summary(mw)
coef(mw, as.survreg = TRUE) ## same interpretation as "survreg"
score_test(mw)
perm_test(mw)
# plot(as.mlt(mw), type = "survivor", newdata = nd1, col = col)



## ----COX-model-fit, echo = FALSE, cache = TRUE--------------------------------
mc <-
Coxph(iDFS ~ randarm, data = CAOsurv, log_first = TRUE)


## ----COX-summary, cache = TRUE, results = "hide", fig.show = "hide"-----------
summary(mc)
score_test(mc)
perm_test(mc)
# plot(as.mlt(mc), type = "survivor", newdata = nd1, col = col)



## ----COX-lHaz, echo = FALSE---------------------------------------------------
## confband
object <- as.mlt(mc)
newdata <- nd1[1,, drop = FALSE]
K <- 20
cheat <- 200
y <- variable.names(object, "response")
q <- mkgrid(object, n = K)[[y]]
q[1] <- q[1] + 1 ## quick fix for log_first
nd <- newdata[rep(1, length(q)), , drop = FALSE]
nd[[y]] <- q
X <- model.matrix(object$model$model, data = nd)
cb <- confint(multcomp::glht(multcomp::parm(coef(object), 
  vcov(object)), linfct = X))$confint
q <- mkgrid(object, n = cheat)[[y]]
q[1] <- q[1] + 1 ## quick fix for log_first 
nd <- newdata[rep(1, length(q)), , drop = FALSE]
nd[[y]] <- q
X <- model.matrix(object$model$model, data = nd)
cb <- confint(multcomp::glht(multcomp::parm(coef(object), 
  vcov(object)), linfct = X), calpha = attr(cb, 
        "calpha"))$confint
cb <- cbind(q, cb)

## ----COX-lHaz-plot------------------------------------------------------------
eval(par_main)
plot(cb[, "q"], cb[, "Estimate"], log = "x", type = "n",
  xlab = lxlab, ylab = ylablHaz, xlim = xlimlHaz <- range(cb[, "q"]),
  ylim = range(cb[, -1]))

polygon(c(cb[, "q"], rev(cb[, "q"])), c(cb[, "lwr"], rev(cb[, "upr"])),
  border = NA, col = rgb(.1, .1, .1, .1))
lines(cb[, "q"], cb[, "Estimate"], lwd = lwd)



## ----STRAT-model-fit, cache = TRUE--------------------------------------------
mcst <- 
Coxph(iDFS | strat ~ randarm, data = CAOsurv, log_first = TRUE)


## ----STRAT-summary, cache = TRUE, results = "hide", fig.show = "hide"---------
summary(mcst)
score_test(mcst)
perm_test(mcst)


## ----STRAT-lHaz-plot----------------------------------------------------------
plot(as.mlt(mcst), newdata = nd2, q = q, type = "trafo", log = "x",
  lty = lty <- 1:4, xlab = lxlab, ylab = ylablHaz, xlim = xlimlHaz,
  col = 1, lwd = lwd)

legend("bottomright", legend = lslevs, title = "Stratum", 
  lty = lty, lwd = lwd, col = 1, bty = "n")




## ----SCOX-model-fit, cache = TRUE---------------------------------------------
mcs <- 
Coxph(iDFS ~ randarm | randarm, data = CAOsurv, log_first = TRUE)


## ----SCOX-summary, cache = TRUE, results = "hide", fig.show = "hide"----------
summary(mcs)
confint(mcs)
perm_test_biv.stram(mcs)
# plot(as.mlt(mcs), type = "survivor", newdata = nd1, col = col)


## ----SCOX-HR-plot, echo = FALSE-----------------------------------------------
qHR <- seq(50, max(q), by = 1)
cumhaz <- predict(mcs, type = "cumhazard", newdata = nd1, q = qHR)
cumhr <- unname(cumhaz[, 2] / cumhaz[, 1])
plot(qHR, cumhr, type = "l", ylab = ylabcumHR, xlab = xlab,
  ylim = ylimHR, xlim = xlimHR <- range(qHR), lwd = lwd)

abline(h = exp(coef(mc)), lty = 2, lwd = 1) ## constant HR
abline(h = 1, lty = 3) ## HR = 1




## ----TCOX-model-fit, cache = TRUE---------------------------------------------
mcv <- 
Coxph(iDFS | randarm ~ 1, data = CAOsurv, log_first = TRUE)

## ----TCOX-summary, cache = TRUE, results = "hide", fig.show = 'hide'----------
logLik(mcv)


## ----TCOX-HR, cache = TRUE----------------------------------------------------
mcv <- as.mlt(mcv)

## grid (was n = 500, qmvnorm failed)
s <- mkgrid(mcv, 40)
s$iDFS <- s$iDFS[s$iDFS >= min(xlimHR) & s$iDFS <= max(xlimHR)]
nd3 <- expand.grid(s)

## confint
K <- model.matrix(mcv$model, data = nd3)
Kyes <- K[nd3$randarm == levels(nd3$randarm)[2],]
Kyes[,grep("Intercept", colnames(Kyes))] <- 0  
gh <- glht(parm(coef(mcv), vcov(mcv)), Kyes)
ci <- exp(confint(gh)$confint)
coxy <- s$iDFS

## confint for constant HR
ci2 <- exp(confint(mc))

## ----TCOX-HR-plot-------------------------------------------------------------
plot(coxy, ci[, "Estimate"], ylim = ylimHR, type = "n",
  xlim = xlimHR, xlab = xlab, ylab = ylabcumHR)
polygon(c(coxy, rev(coxy)), c(ci[,"lwr"], rev(ci[, "upr"])),
        border = NA, col = rgb(.1, .1, .1, .1))
lines(coxy, ci[, "Estimate"], lty = 1, lwd = lwd)

## constant HR
polygon(c(coxy[c(1, length(coxy))], rev(coxy[c(1, length(coxy))])),
  rep(ci2, c(2, 2)), border = NA, col = rgb(.1, .1, .1, .1))
abline(h = exp(coef(mc)), lty = 2, lwd = 1)

## HR = 1
abline(h = 1, lty = 3)


## ----DEPCENS-preproc, echo = FALSE--------------------------------------------
## DepC: loss of follow-up (everyone else is admin censored) Mail TH 23-06-12
patnr_lofu <-c(1012, 2003, 3002, 3003, 6018, 7001, 7003, 7005, 7008, 7012, 10003,
              10012, 11018, 12003, 12014, 13028, 14002, 15001, 16001, 16004, 16005,
              16007, 16009, 18016, 18025, 21011, 21013, 21014, 21022, 21023, 21026,
              21027, 21029, 21043, 22003, 23008, 24008, 24021, 25001, 25004, 25005,
              25006, 26005, 26018, 27005, 27030, 27034, 29002, 30006, 30011, 31003,
              31004, 31005, 34001, 35011, 35014, 36017, 41004, 42001, 42003, 42005,
              42007, 42010, 44004, 44005, 45002, 45003, 45009, 45011, 46003, 49001,
              49003, 49011, 49012, 49015, 50001, 50003, 50004, 50007, 50011, 52004,
              54004, 56006, 56008, 59002, 59005, 68001, 70010, 71002, 73009, 74004,
              75002, 75004, 75005, 80003, 81001, 84005, 84007, 86002) 
ilofu <- with(CAOsurv, which(patnr %in% patnr_lofu))
CAOsurv$DepCevent <- CAOsurv$OSevent
CAOsurv$DepCevent <- factor(as.numeric(CAOsurv$DepCevent), levels = 0:2,
  labels = c("AdminC", "EoI", "DepC"))
CAOsurv$DepCevent[ilofu] <- "DepC"


## ----DEPCENS-table, results = 'hide'------------------------------------------
CAOsurv$nDepCevent <- factor(as.character(CAOsurv$DepCevent),
  levels = c("AdminC", "EoI", "DepC"), 
  labels = c("Administrative censoring", "Event of interest", "Loss of follow-up"))
tab <- xtabs(~ nDepCevent + randarm, data = CAOsurv)
tab




## ----DEPCENS-model-fit, cache = TRUE------------------------------------------
md <- 
Coxph(Surv(OStime, event = DepCevent) ~ randarm, data = CAOsurv)


## ----DEPCENS-summary, cache = TRUE, results = "hide", fig.show = 'hide'-------
summary(md)
confint(md)


## ----COXME-install, echo = TRUE, eval = FALSE---------------------------------
## install.packages("tramME")
## library("tramME")


## ----COXME-load---------------------------------------------------------------
library("tramME")



## ----COXME-model-fit, cache = TRUE--------------------------------------------
mcME <- 
CoxphME(iDFS ~ randarm + (1 | Block), data = CAOsurv,
  log_first = TRUE)


## ----COXME-summary, cache = TRUE, results = "hide", fig.show = 'hide'---------
summary(mcME)
confint(mcME)




## ----COXME-margsurv, eval = FALSE---------------------------------------------
## ## computationally intensive
## if (!file.exists("ME-margdist.rda")) {
## mod <- mcME
## 
## ## A function to evaluate the joint cdf of the response and the random effects:
## ## Takes a vector of random effect and covariates values, evaluates the conditional
## ## distribution at these values and multiplies it with the pdf of the random effects
## jointCDF <- function(re, nd, mod) {
## nd <- nd[rep(1, length(re)), ]
## nd$Block <- seq(nrow(nd)) ## to take vector-valued REs
## pr <- predict(mod, newdata = nd, ranef = re, type = "distribution") *
## dnorm(re, 0, sd = sqrt(varcov(mod)[[1]][1, 1]))
## c(pr)
## }
## ## Marginalize the joint cdf by integrating out the random effects
## ## using adaptive quadrature
## marginalCDF <- function(nd, mod) {
## nd$cdf <- integrate(jointCDF, lower = -Inf, upper = Inf, nd = nd, mod = mod)$value
## nd
## }
## ## Set up the grid on which we evaluate the marginal distribution
## nd <- expand.grid(iDFS = 1:max(CAOsurv$DFStime), randarm = unique(CAOsurv$randarm))
## ## Calls marginalCDF on each row of nd
## ## (done in parallel to speed up computations)
## mp <- parallel::mclapply(split(nd, seq(nrow(nd))),
##   marginalCDF, mod = mod, mc.cores = 4)
## mp <- do.call("rbind", mp)
## save(mp, file = "ME-margdist.rda")
## } else load("ME-margdist.rda")
## mp$surv <- with(mp, 1 - cdf)
## 
## plot(surv_iDFS, ylim = ylimS, xlim = xlim,
##   col = lcol, lwd = lwd, xlab = xlab, ylab = ylabS)
## with(mp[mp$randarm == levs[1], ], lines(iDFS, surv, col = col[1], lwd = lwd))
## with(mp[mp$randarm == levs[2], ], lines(iDFS, surv, col = col[2], lwd = lwd))
## legend("bottomright", legend = levs, col = col, bty = "n", lty = 1, lwd = 1, cex = .8)


## ----MCOX-preproc, echo = FALSE-----------------------------------------------
### convert "exact" event dates to interval-censoring (+/- two days)
tmp <- CAOsurv$iDFS
exact <- tmp[, 3] == 1
tmp[exact, 2] <- tmp[exact, 1] + 2
tmp[exact, 1] <- pmax(tmp[exact, 1] - 2, 0)
tmp[exact, 3] <- 3
CAOsurv$iDFS2 <- tmp


## ----MCOX-model-fit, cache = TRUE---------------------------------------------
mmc <- 
mtram(Coxph(iDFS2 ~ randarm, data = CAOsurv, log_first = TRUE),
  formula = ~ (1 | Block), data = CAOsurv)


## ----MCOX-FUN-----------------------------------------------------------------
## marginal HR from "mtram"
## <FIXME> reset seed on.exit </FIXME>
mHR.mtram <- function(object, with_confint = FALSE, seed = 1) {
  stopifnot(inherits(object, "mtram"))
  cf <- coef(object)
  cf <- cf[-grep("Bs", names(cf))]
  stopifnot(length(cf) == 2)
  mlHR <- cf[1] / sqrt(1 + cf["gamma1"]^2)
  ret <- mHR <- exp(mlHR)
  if (with_confint) {
    set.seed(seed)
    S <- vcov(object)
    rbeta <- rmvnorm(10000, mean = coef(object), sigma = S)
    s <- rbeta[,ncol(rbeta)]
    rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
    ci <- quantile(exp(rbeta[, ncol(rbeta)]), prob = c(.025, .975))
    ret <- c(mHR, ci)
    ret <- as.array(t(ret))
  }
  return(ret)
}

## ----MCOX-summary, cache = TRUE, results = "hide", fig.show = 'hide'----------
coef(mmc)
sqrt(diag(vcov(mmc)))
(ci_MCOX <- mHR.mtram(mmc, with_confint = TRUE))



## ----HTECOX-model-fit, cache = TRUE-------------------------------------------
ma <- 
CoxphME(iDFS ~ randarm + s(age, by = as.ordered(randarm),
    fx = TRUE, k = 6), data = CAOsurv, log_first = TRUE)
nd <- model.frame(ma)[rep(2, 100), ]
nd$age <- seq(min(CAOsurv$age), max(CAOsurv$age), length.out = 100)
xx <- model.matrix(ma, data = nd, type = "X", keep_sign = FALSE)$X
ip <- grep("randarm", names(bb <- coef(ma, with_baseline = TRUE)))
vc <- vcov(ma, parm = names(bb)[ip])
bb <- bb[ip]

## NOTE: unadjusted
cb <- exp(confint(multcomp::glht(multcomp::parm(bb, vc), linfct = xx),
                  calpha = univariate_calpha())$confint)

## ----HTECOX-summary, cache = TRUE, results = "hide", fig.show = 'hide'--------
summary(ma)


## ----HTECOX-HR-plot-----------------------------------------------------------
## Plot HR
plot(nd$age, cb[, "Estimate"], type = "n", ylab = "Hazard ratio", xlab = "Age (in years)",
     ylim = ylimHR)
polygon(c(nd$age, rev(nd$age)), c(cb[, "lwr"], rev(cb[, "upr"])),
        border = NA, col = rgb(.1, .1, .1, .1))
lines(nd$age, cb[, "Estimate"], lwd = lwd)
abline(h = 1, lty = 3)
rug(CAOsurv$age, lwd = 2, col = rgb(.1, .1, .1, .1))


## ----TRT-load-install, echo = TRUE, eval = FALSE------------------------------
## install.packages("trtf")
## library("trtf")

## ----TRT-load-----------------------------------------------------------------
library("trtf")
set.seed(4)


## ----TRTF-model-fit, cache = TRUE, warning=TRUE-------------------------------
tr <-
trafotree(Coxph(iDFS ~ randarm, data = CAOsurv, log_first = TRUE),
  formula = iDFS ~ randarm | age, data = CAOsurv,
  control = ctree_control(teststat = "maximum", minbucket = 40))

## ----TRT-results, results = "hide", fig.show = 'hide'-------------------------
logLik(tr)


## ----TRTF-surv-plot, fig.width = 10, fig.height = 6---------------------------
library("ATR")
plot(rotate(tr), tp_args = list(newdata = nd1, type = "survivor", col = col, lwd = lwd),
  terminal_panel = trtf:::node_mlt)



## ----FRAILTY-model-fit, cache = TRUE------------------------------------------
mf <- 
Coxph(iDFS ~ randarm, data = CAOsurv, log_first = TRUE, frailty = "Gamma")


## ----FRAILTY-summary, cache = TRUE, results = "hide", fig.show = 'hide'-------
logLik(mf)
coef(mf)[trt]
coef(mf, addparm = TRUE)
confint(mf, parm = c(trt, "logrho"))



## ----LOGIT-model-fit----------------------------------------------------------
ml <- 
Colr(iDFS ~ randarm, data = CAOsurv, log_first = TRUE)


### Appendix
## ----pkgs---------------------------------------------------------------------
## additional packages
pkgs <- c("fastGHQuad", "icenReg", "TransModel", "rms", "ICsurv", "eha",
  "rstpm2", "flexsurv", "mpr", "gamlss", "gamlss.cens", 
  "coxme", "parfm", "frailtyEM", "frailtypack", "mgcv", "timereg")

## ----install-pkgs-------------------------------------------------------------
ix <- which(!sapply(pkgs, require, char = TRUE))
if (length(ix) > 0) {install.packages(pkgs[ix], repos = "https://stat.ethz.ch/CRAN/")
   sapply(pkgs[ix], require, char = TRUE)}


## ----pkgs-setup---------------------------------------------------------------
`coef<-` <- mlt::`coef<-` ## masked by pkg rstpm2
Surv <- survival::Surv ## masked by eha



## ----result-summary, include = FALSE------------------------------------------
frmtcall <- function(mod, call = NA, tex = TRUE) {
  pkg <- NA
  ## if call is specified
  if (!is.na(call)) {
  call <- rev(strsplit(call, "::")[[1]])
  pkg <- call[2]
  call <- call[1]
  } else {
  if (inherits(mod, "emfrail")) call <- attr(mod, "call")[[1]]
  else call <- mod$call[[1]]
  call <- as.character(call)
  if (length(call) > 1) {
    pkg <- call[2]
    call <- call[3]
  }}
  fct <- ifelse(tex, paste0("\\code{", gsub("_", "\\\\_", call), "}"), call)
  ret <- c("Call" = call, "Function" = fct)
  if (!is.na(pkg)) {
    pkg <- ifelse(tex, paste0("\\pkg{", gsub("_", "\\\\_", pkg), "}"), pkg)
    ret <- c("Call" = call, "Function" = fct, "Package" = pkg)
  }
  return(ret)
}

logLik.aftreg <- logLik.phreg <- function(object) object$loglik[2]
logLik.stpm2 <- function(object) -object@minuslogl(object@coef)
logLik.mpr <- function(object) object$model$loglike
logLik.icenReg_fit <- function(object) summary(object)$llk
logLik.cox.aalen <- function(object) NA ## no reported logLik found
vcov.aftreg <- vcov.phreg <- function(object) object$var
vcov.stpm2 <- function(object) object@vcov
vcov.emfrail <- function(object) {
  vc <- object$var
  ncf <- names(object$coefficients)
  colnames(vc) <- c(ncf, (length(ncf) + 1):ncol(vc))
  rownames(vc) <- colnames(vc)
  return(vc)
}
logLik.frailtyPenal <- function(object) object$logLikPenal ## pen. marg. logLik
vcov.frailtyPenal <- function(object) {
  vc <- matrix(object$varH)
  dimnames(vc) <- list(names(object$coef), names(object$coef))
  return(vc)
}
coef.frailtyPenal <- function(object) object$coef

link.stpm2.po <- function (S) -logit(as.vector(S))
link.stpm2.ph <- function (S) log(-log(as.vector(S)))

tab <- function(mod, parm = trt, math = TRUE, mark = TRUE, tex = TRUE) {
  ll <- logLik(mod)
  
  if (inherits(mod, "mpr")) cf <- coef(mod)[[1]][parm <- paste(parm, "b", sep = ".")]
  else cf <- coef(mod)
  if (inherits(mod, "cox.aalen")) cf <- cf[, "Coef."]
  cf <- ifelse(length(cf) > 1, cf[parm], cf)
  if (inherits(mod, "tram")) cf <- c(1, -1)[mod$negative + 1] * cf
  
  if (inherits(mod, "ic_ph") | inherits(mod, "ic_po")) se <- ifelse(math, "$-$", NA)
  else {
    vc <- vcov(mod)
    se <- ifelse(sum(dim(vc)) > 2, sqrt(vc[parm, parm]), sqrt(vc))
  }
  
  call <- NA
  if (inherits(mod, "stpm2")) call <- "rstpm2::stpm2"
  if (inherits(mod, "cox.aalen")) call <- "timereg::Gprop.odds"
  call <- frmtcall(mod, call = call, tex = tex)
  
  cfint <- switch(call["Call"],
    "coxph" = "log-HR", "Coxph" = "log-HR", "CoxphME" = "log-HR", 
    "gam" = ifelse(mod$family$family == "Cox PH", "log-HR", NA),
    "coxme" = "log-HR", "frailtyPenal" = "log-HR", "emfrail" = "log-HR",
    "survreg" = "log-AF", "Survreg" = "log-HR",
    "stpm2" = ifelse(identical(mod@link$link, link.stpm2.ph, ignore.environment = TRUE),
      "log-HR", ifelse(identical(mod@link$link, link.stpm2.po, ignore.environment = TRUE), 
        "log-OR", NA)),
    "flexsurvspline" = switch(mod$scale, "hazard" = "log-HR", "odds" = "$-$log-OR"), 
    "flexsurvreg" = switch(mod$call$dist, "weibullPH" = "log-HR", "weibull" = "log-AF"),
    "ic_par" = switch(mod$call$model, "ph" = "log-HR", "aft" = "log-AF"), 
    "ic_sp" = switch(mod$call$model, "ph" = "log-HR", "po" = "$-$log-OR", "aft" = "log-AF"), 
    "cph" = "log-HR", "phreg" = "log-HR", "Gprop.odds" = "log-OR", "Colr" = "log-OR")
  if (inherits(mod, c("stram", "mpr", "gamlss"))) cfint <- NA
  if (is.null(cfint)) cfint <- NA
  
  mark <- ifelse(mark && call["Call"] %in% c("gam", "coxph", "coxme", "cph", "ic_sp", "ic_ph", "ic_po", 
    "cox.aalen", "emfrail", "frailtyPenal"), TRUE, FALSE)
  ret <- c(cfint, frmt3(cf, math = math), frmt3(se, math = math),
    frmtll(ll, mark = mark, math = math))
  names(ret) <- n <- c("Interpretation", "Estimate", "Std. Error", "logLik")
  c(call[-1], ret)
}

print.results <- function(objects) {
  ret <- lapply(objects, function(x) tab(x, math = FALSE, mark = FALSE, tex = FALSE))
  data.frame(do.call("rbind", ret), check.names = FALSE)
}




## ----WEI-iDFS-fit, cache = TRUE-----------------------------------------------
mwi1 <- tram::Survreg(iDFS ~ randarm, data = CAOsurv, dist = "weibull")
mwi2 <- icenReg::ic_par(iDFS ~ randarm, data = CAOsurv, dist = "weibull",
  model = "ph")
mwi3 <- flexsurv::flexsurvreg(iDFS ~ randarm, data = CAOsurv,
  dist = "weibullPH")
mwi4 <- survival::survreg(iDFS ~ randarm, data = CAOsurv, dist = "weibull")
mwi5 <- icenReg::ic_par(iDFS ~ randarm, data = CAOsurv, dist = "weibull",
  model = "aft") 


## ----WEI-iDFS-results---------------------------------------------------------
print.results(list(mwi1, mwi2, mwi3, mwi4, mwi5))


## ----preproc-interval, results = "hide"---------------------------------------
## is needed for "rstpm2" package
CAOsurv$iDFSevent <- as.numeric(CAOsurv$DFSevent)
ic <- with(CAOsurv, which(is.finite(iDFStime2) &  iDFStime2 > iDFStime))
CAOsurv$iDFSevent[ic] <- 3
table(CAOsurv$iDFSevent)
with(CAOsurv, all.equal(Surv(time = iDFStime, time2 = iDFStime2, event = iDFSevent,
    type = "interval"), iDFS)) ## check


## ----Cox-iDFS-fit, echo = FALSE, cache = TRUE---------------------------------
mci1 <- tram::Coxph(iDFS ~ randarm, data = CAOsurv, log_first = TRUE)
mci2 <- rstpm2::stpm2(Surv(time = iDFStime, time2 = iDFStime2,
  event = iDFSevent, type = "interval") ~ randarm, data = CAOsurv)
mci3 <- flexsurv::flexsurvspline(iDFS ~ randarm, data = CAOsurv, k = 3)
mci4 <- icenReg::ic_sp(iDFS ~ randarm, data = CAOsurv, model = "ph")


## ----Cox-iDFS-tab-print-------------------------------------------------------
print.results(list(mci1, mci2, mci3, mci4))




## ----Cox-DFS-fit, echo = FALSE, cache = TRUE----------------------------------
mc1 <- tram::Coxph(DFS ~ randarm, data = CAOsurv, log_first = TRUE)
mc2 <- survival::coxph(DFS ~ randarm, data = CAOsurv)
mc3 <- rms::cph(DFS ~ randarm, data = CAOsurv)


## ----Cox-DFS-results, eval = TRUE---------------------------------------------
print.results(list(mc1, mc2, mc3))








## ----STRAT-iDFS-fit, cache = TRUE---------------------------------------------
mstci1 <- tram::Coxph(iDFS | strat ~ randarm, data = CAOsurv, log_first = TRUE)
mstci2 <- rstpm2::stpm2(Surv(time = iDFStime, time2 = iDFStime2, event = iDFSevent,
    type = "interval") ~ randarm + strata(strat), data = CAOsurv)
mstci3 <- NULL
#mstci3 <- flexsurv::flexsurvspline(iDFS ~ randarm + gamma1(strat) + gamma2(strat),
#  data = CAOsurv, k = 3, control=list(ndeps=rep(1e-09,12)))
## COMMENT: "flexsurv" package: look at plot(model)
## control: necessary for flexsurv 2.3; see email by Nov 4
##          by Chris Jackson
##          however, the hessian is singular, so we exclude this model
##          for the time being.


## ----STRAT-iDFS-results-------------------------------------------------------
print.results(list(mstci1, mstci2))##, mstci3))



## ----STRAT-DFS-fit, cache = TRUE----------------------------------------------
mstc1 <- tram::Coxph(DFS | strat ~ randarm, data = CAOsurv, log_first = TRUE)
mstc2 <- survival::coxph(DFS ~ randarm + strata(strat), data = CAOsurv)
mstc3 <- rms::cph(DFS ~ randarm + strat(strat), data = CAOsurv)


## ----STRAT-DFS-results--------------------------------------------------------
print.results(list(mstc1, mstc2, mstc3))



## ----STRAT-Wei-iDFS-fit, echo = FALSE, cache = TRUE---------------------------
mstw1 <- tram::Survreg(DFS | strat ~ randarm, data = CAOsurv)
mstw2 <- eha::phreg(DFS ~ randarm + strata(strat), data = CAOsurv)
mstw3 <- survival::survreg(DFS ~ randarm + strata(strat), data = CAOsurv)


## ----STRAT-Wei-iDFS-results---------------------------------------------------
print.results(list(mstw1, mstw2, mstw3))




## ----LS-iDFS-Wei-fit, cache = TRUE--------------------------------------------
mswi1 <- tram::Survreg(iDFS ~ randarm | randarm, data = CAOsurv,
  remove_intercept = FALSE)
tmp <- CAOsurv[, c("iDFS", "randarm")] ## NA in other columns prompts error
gen.cens(family = "WEI2", type = "interval")
mswi2 <- gamlss::gamlss(formula = iDFS ~ randarm, sigma.fo = ~ randarm,
  family = gamlss.cens::cens(family = "WEI2", type = "interval"),
  data = tmp, control = gamlss.control(n.cyc = 300, trace = FALSE))



## ----LS-iDFS-results----------------------------------------------------------
print.results(list(mswi1, mswi2))




## ----LS-DFS-fit, echo = FALSE, cache = TRUE-----------------------------------
msw1 <- tram::Survreg(DFS ~ randarm | randarm, data = CAOsurv,
  remove_intercept = FALSE)
msw2 <- mpr::mpr(DFS ~ list(~ randarm, ~ randarm), data = CAOsurv)



## ----LS-DFS-results-----------------------------------------------------------
print.results(list(msw1, msw2))




## ----TVAR-iDFS-fit, cache = TRUE----------------------------------------------
mcvi1 <- tram::Coxph(iDFS | randarm ~ 1, data = CAOsurv)
mcvi2 <- flexsurv::flexsurvspline(iDFS ~ randarm +
    gamma1(randarm) + gamma2(randarm), data = CAOsurv, k = 3)



## ----TVAR-iDFS-plot, fig.width = 6, fig.height = 3----------------------------
## cumHR from "tram"
xlim.tvar <- c(100, max(q))

y <- variable.names(mcvi1, "response")
s <- mkgrid(mcvi1, n = 50)
s[[y]] <- s[[y]][s[[y]] > xlim[1] & s[[y]] < xlim[2]]

cumhaz <- predict(as.mlt(mcvi1), newdata = s, type = "cumhazard")
cumhr <- cumhaz[,2] / cumhaz[,1]

par(mgp = c(2.5, 1, 0), mar = c(4, 4, 1.5, 4))
plot(s[[y]], cumhr, ylim = ylimHR, type = "l",
     xlab = xlab, ylab = ylabcumHR, las = 1, lwd = lwd)
abline(h = 1, lty = 3)

## cumHR from "flexsurvspline"
cumhaz <- predict(mcvi2, type = "cumhaz", newdata =  nd1)
cumhr <- unlist(unname(cumhaz[[1]][[2]][2] /  cumhaz[[1]][[1]][2]))
t <- unlist(unname(cumhaz[[1]][[1]][1]))
ret <- as.data.frame(cbind(t, cumhr))
ret <- ret[ret$t > xlim.tvar[1] & ret$t < xlim.tvar[2], ]
lines(ret$t, ret$cumhr, lty = 2, lwd = 2, col = col2 <- "darkgrey")
legend("topright", lty = 1:2, lwd = 2, col = c("black", col2),
  legend = c(bquote("package:"~bold("tram")), bquote("package:"~bold("flexsurv"))),
  bty = "n")



## ----TVAR-DFS-fit, cache = TRUE-----------------------------------------------
mcv1 <- tram::Coxph(DFS | randarm ~ 1, data = CAOsurv, log_first = TRUE)
mcv2 <- flexsurv::flexsurvspline(DFS ~ randarm + gamma1(randarm) + gamma2(randarm),
  data = CAOsurv, k = 3)


## ----TVAR-DFS-plot, fig.width = 6, fig.height = 3-----------------------------
## cumHR from "tram"
xlim.tvar <- c(100, max(q))

y <- variable.names(mcv1, "response")
s <- mkgrid(mcv1, n = 50)
s[[y]] <- s[[y]][s[[y]] > xlim[1] & s[[y]] < xlim[2]]

cumhaz <- predict(as.mlt(mcv1), newdata = s, type = "cumhazard")
cumhr <- cumhaz[,2] / cumhaz[,1]

par(mgp = c(2.5, 1, 0), mar = c(4, 4, 1.5, 4))
plot(s[[y]], cumhr, ylim = ylimHR, type = "l",
     xlab = xlab, ylab =  ylabcumHR, las = 1, lwd = lwd)
abline(h = 1, lty = 3)

## cumHR from "flexsurvspline"
cumhaz <- predict(mcv2, type = "cumhaz", newdata =  nd1)
cumhr <- unlist(unname(cumhaz[[1]][[2]][2] /  cumhaz[[1]][[1]][2]))
t <- unlist(unname(cumhaz[[1]][[1]][1]))
ret <- as.data.frame(cbind(t, cumhr))
ret <- ret[ret$t > xlim.tvar[1] & ret$t < xlim.tvar[2], ]
lines(ret$t, ret$cumhr, lty = 2, lwd = 2, col = col2 <- "darkgrey")
legend("topright", lty = 1:2, lwd = 2, col = c("black", col2),
  legend = c(bquote("package:"~bold("tram")), bquote("package:"~bold("flexsurv"))),
  bty = "n")



## ----MIXED-DFS-fit, cache = TRUE----------------------------------------------
mcME1 <- tramME::CoxphME(DFS ~ randarm + (1 | Block), data = CAOsurv, log_first = TRUE)
mcME2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv,
  cluster = "Block", RandDist = "LogN")
mcME3 <- coxme::coxme(DFS ~ randarm + (1 | Block), data = CAOsurv)


## ----MIXED-DFS-results--------------------------------------------------------
print.results(list(mcME1, mcME2, mcME3))




## ----HTECOX-DFS-fit-----------------------------------------------------------
ma1 <- CoxphME(DFS ~ randarm +
    s(age, by = as.ordered(randarm), fx = TRUE, k = 6),
               data = CAOsurv, log_first = TRUE)
ma2 <- gam(DFStime ~ randarm +
    s(age, by = as.ordered(randarm), fx = TRUE, k = 6),
               data = CAOsurv, family = cox.ph(), weights = DFSevent)


## ----HTECOX-DFS-results-------------------------------------------------------
print.results(list(ma1, ma2))


## ----HTECOX-DFS-plot----------------------------------------------------------
nd <- model.frame(ma1)[rep(2, 100), ]
nd$age <- seq(min(CAOsurv$age), max(CAOsurv$age), length.out = 100)
xx <- model.matrix(ma1, data = nd, type = "X", keep_sign = FALSE)$X
ip <- grep("randarm", names(bb <- coef(ma1, with_baseline = TRUE)))
vc <- vcov(ma1, parm = names(bb)[ip])
bb <- bb[ip]

cb1 <- exp(confint(multcomp::glht(multcomp::parm(bb, vc), linfct = xx),
                  calpha = univariate_calpha())$confint)

plot(nd$age, cb1[, "Estimate"], type = "n", ylab = "Hazard ratio", xlab = "Age (in years)",
     ylim = ylimHR)
matlines(nd$age, cb1, lwd = lwd, col = 1, lty = 1)
abline(h = 1, lty = 3)

summary(ma2)

nd2 <- model.frame(ma2)[rep(2, 100), ]
nd2$age <- seq(min(CAOsurv$age), max(CAOsurv$age), length.out = 100)
pr <- predict(ma2, newdata = nd2, type = "link", se.fit = TRUE)

matlines(nd2$age, exp(c(pr$fit) + qnorm(0.975) * pr$se.fit %o% c(0, -1, 1)),
         col = col2, lwd = lwd, lty = 2)

legend("bottomright", lty = 1:2, lwd = 2, col = c("black", col2),
  legend = c(bquote("package:"~bold("tramME")), bquote("package:"~bold("mgcv"))),
  bty = "n")



## ----FRAILTY-DFS-fit, cache = TRUE--------------------------------------------
mfc1 <- tram::Coxph(DFS ~ randarm, data = CAOsurv, log_first = TRUE, frailty = "Gamma")
mfc2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv,
  cluster = "id", RandDist = "Gamma")
mfc3 <- survival::coxph(DFS ~ randarm + frailty(id, distribution = "gamma"), data = CAOsurv)
mfc4 <- frailtyEM::emfrail(DFS ~ randarm + cluster(id), data = CAOsurv)
mfc5 <- frailtypack::frailtyPenal(DFS ~ randarm + cluster(id), data = CAOsurv,
  RandDist = "Gamma", n.knots = 10, kappa = 1)


## ----FRAILTY-DFS-results------------------------------------------------------
print.results(list(mfc1, mfc2, mfc3, mfc4, mfc5))




## ----Colr-DFS-fit, cache = TRUE-----------------------------------------------
mo1 <- tram::Colr(DFS ~ randarm, data = CAOsurv, log_first = TRUE)
mo2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv, link.type = "PO")
mo3 <- flexsurv::flexsurvspline(iDFS ~ randarm, data = CAOsurv, k = 3, scale = "odds")
mo4 <- timereg::Gprop.odds(DFS ~ prop(randarm), data = CAOsurv)


## ----Colr-DFS-results---------------------------------------------------------
print.results(list(mo1, mo2, mo3, mo4))


## ----session, results = "markup"----------------------------------------------
sessionInfo()

