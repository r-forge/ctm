### Code from
### "Smooth Transformation Models for Survival Analysis: A Tutorial Using R"
###   by Sandra Siegfried, Balint Tamasi & Torsten Hothorn

library("tram")

x <- vignette("survtram", package = "tram")
source(file.path(x$Dir, "doc", x$R), echo = TRUE)

### Appendix

## ----pkgs---------------------------------------------------------------------
## additional packages
pkgs <- c("fastGHQuad", "icenReg", "TransModel", "rms", "ICsurv", "eha",
  "rstpm2", "flexsurv", "mpr", "gamlss", "gamlss.cens", 
  "coxme", "parfm", "frailtyEM", "frailtypack", "mgcv", "timereg")

## ----install-pkgs-------------------------------------------------------------
ix <- which(!sapply(pkgs, require, char = TRUE))
if (length(ix) > 0) {
   install.packages(pkgs[ix], repos = "https://stat.ethz.ch/CRAN/")
   sapply(pkgs[ix], require, char = TRUE)
}


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
mci1 <- tram::Coxph(iDFS ~ randarm, data = CAOsurv)
mci2 <- rstpm2::stpm2(Surv(time = iDFStime, time2 = iDFStime2,
  event = iDFSevent, type = "interval") ~ randarm, data = CAOsurv)
mci3 <- flexsurv::flexsurvspline(iDFS ~ randarm, data = CAOsurv, k = 3)
mci4 <- icenReg::ic_sp(iDFS ~ randarm, data = CAOsurv, model = "ph")


## ----Cox-iDFS-tab-print-------------------------------------------------------
print.results(list(mci1, mci2, mci3, mci4))




## ----Cox-DFS-fit, echo = FALSE, cache = TRUE----------------------------------
mc1 <- tram::Coxph(DFS ~ randarm, data = CAOsurv)
mc2 <- survival::coxph(DFS ~ randarm, data = CAOsurv)
mc3 <- rms::cph(DFS ~ randarm, data = CAOsurv)


## ----Cox-DFS-results, eval = TRUE---------------------------------------------
print.results(list(mc1, mc2, mc3))








## ----STRAT-iDFS-fit, cache = TRUE---------------------------------------------
mstci1 <- tram::Coxph(iDFS | strat ~ randarm, data = CAOsurv)
mstci2 <- rstpm2::stpm2(Surv(time = iDFStime, time2 = iDFStime2, event = iDFSevent,
    type = "interval") ~ randarm + strata(strat), data = CAOsurv)
# mstci3 <- flexsurv::flexsurvspline(iDFS ~ randarm + gamma1(strat) + gamma2(strat),
#  data = CAOsurv, k = 3, control=list(ndeps=rep(1e-09,12)))
## control: necessary for flexsurv 2.3; see email by Nov 4
##          by Chris Jackson
##          however, the hessian is singular, so we exclude this model
##          for the time being.


## ----STRAT-iDFS-results-------------------------------------------------------
print.results(list(mstci1, mstci2))##, mstci3))



## ----STRAT-DFS-fit, cache = TRUE----------------------------------------------
mstc1 <- tram::Coxph(DFS | strat ~ randarm, data = CAOsurv)
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
mcv1 <- tram::Coxph(DFS | randarm ~ 1, data = CAOsurv)
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
mcME1 <- tramME::CoxphME(DFS ~ randarm + (1 | Block), data = CAOsurv)
mcME2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv,
  cluster = "Block", RandDist = "LogN")
mcME3 <- coxme::coxme(DFS ~ randarm + (1 | Block), data = CAOsurv)


## ----MIXED-DFS-results--------------------------------------------------------
print.results(list(mcME1, mcME2, mcME3))




## ----HTECOX-DFS-fit-----------------------------------------------------------
ma1 <- CoxphME(DFS ~ randarm +
    s(age, by = as.ordered(randarm), fx = TRUE, k = 6),
               data = CAOsurv)
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
mfc1 <- tram::Coxph(DFS ~ randarm, data = CAOsurv, frailty = "Gamma")
mfc2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv,
  cluster = "id", RandDist = "Gamma")
mfc3 <- survival::coxph(DFS ~ randarm + frailty(id, distribution = "gamma"), data = CAOsurv)
mfc4 <- frailtyEM::emfrail(DFS ~ randarm + cluster(id), data = CAOsurv)
mfc5 <- frailtypack::frailtyPenal(DFS ~ randarm + cluster(id), data = CAOsurv,
  RandDist = "Gamma", n.knots = 10, kappa = 1)


## ----FRAILTY-DFS-results------------------------------------------------------
print.results(list(mfc1, mfc2, mfc3, mfc4, mfc5))




## ----Colr-DFS-fit, cache = TRUE-----------------------------------------------
mo1 <- tram::Colr(DFS ~ randarm, data = CAOsurv)
mo2 <- rstpm2::stpm2(Surv(DFStime, DFSevent) ~ randarm, data = CAOsurv, link.type = "PO")
mo3 <- flexsurv::flexsurvspline(iDFS ~ randarm, data = CAOsurv, k = 3, scale = "odds")
mo4 <- timereg::Gprop.odds(DFS ~ prop(randarm), data = CAOsurv)


## ----Colr-DFS-results---------------------------------------------------------
print.results(list(mo1, mo2, mo3, mo4))


## ----session, results = "asis"------------------------------------------------
toLatex(sessionInfo(), locale = FALSE)

