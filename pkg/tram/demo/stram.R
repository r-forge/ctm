### Demo for location-scale transformation models ###
### DOI: 10.1080/00031305.2023.2203177

## run all the code (including the computationally intensive examples)
run_all <- FALSE 

## general setup ##
library("lattice")
library("gridExtra")
library("grid") ### for textGrob
library("latticeExtra")
library("reshape2")
library("colorspace")
library("xtable")

acol <- sequential_hcl(6, "BluYl")[1:5]

col <- acol[c(2, (length(acol)) - 1)]
lcol <- lighten(col, amount = .4) ## lighten color for overlayed lines 

## lattice
trellis.par.set(
  list(
    plot.symbol = list(col = 1, pch = 20, cex = 0.7),
    box.rectangle = list(col = 1),
    box.umbrella = list(lty = 1, col = 1),
    strip.background = list(col = "white")
  )
)

ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

## Formatting
big.mark <- ","
frmt0 <- round

is.neg <- function(x) x < 0

frmt1 <- function(x, math = TRUE) {
    if (!is.numeric(x)) return(x)
    ret <- formatC(round(x, 1), digits = 1, format = "f", big.mark = big.mark) 
    if (math) ret <- paste("$", ret, "$")
    if (is.matrix(x)) {
        ret <- matrix(ret, nrow = nrow(x))
        dimnames(ret) <- dimnames(x)
    }
    ret
}

frmt4 <- function(x, math = TRUE) {
    if (!is.numeric(x)) return(x)
    ret <- formatC(round(x, 4), digits = 4, format = "f", big.mark = big.mark) 
    if (math) ret <- paste("$", ret, "$")
    if (is.matrix(x)) {
        ret <- matrix(ret, nrow = nrow(x))
        dimnames(ret) <- dimnames(x)
    }
    ret
}

frmtI <- function(x, math = TRUE) {
  if (is.character(x)) return(x)
  ret <- trimws(formatC(x, format = "fg", width = 7, big.mark = big.mark))
  if (x < 0) ret <- paste0(ifelse(math, "$-$", "-"), ret)
  if (math) ret <- paste("$", ret, "$")
  if (is.matrix(x)) {
    ret <- matrix(ret, nrow = nrow(x))
    dimnames(ret) <- dimnames(x)
  }
  ret
}
    
toUpper <- function(x) gsub("\\b([[:lower:]])([[:lower:]]+)", "\\U\\1\\L\\2",
                            x, perl = TRUE)

frmtll <- function(ll) paste("logLik = ", frmt1(ll, math = FALSE))

## definitions
tram <- "Location transformation model"
stram <- "Location-scale transformation model"


## Section 3.1.1. Stratification ##
## ----STRAT-setup, include = FALSE---------------------------------------------
## libraries
library("tram")
library("survival")

## ----STRAT-preproc, include = FALSE-------------------------------------------
### blood loss data: DOI: 10.1111/jth.14795 
load(system.file("rda", "bloodloss.rda", package = "TH.data"))

## delivery mode
blood$mode <- with(blood, SECTIO.prim == "yes" | SECTIO.sek == "yes" | SECTIO.not == "yes")
blood$mode <- factor(blood$mode, levels = c(FALSE, TRUE),
                     labels = c("Vaginal delivery", "Cesarean section"))

## interval censoring for appropriate likelihood
sMBL <- sort(unique(blood$MBL))
blood$MBLc <- cut(blood$MBL, breaks = c(-Inf, sMBL), ordered = TRUE)
i <- unclass(blood$MBLc)
blood$MBLi <- Surv(c(0, sMBL[-length(sMBL)])[i], sMBL[i], type = "interval2")

## remove outlier
blood <- blood[blood$MBL < 4000, ]
(tab <- table(blood$mode))

## ----STRAT-model--------------------------------------------------------------
OR <- 15

## stratified transformation model
m <- BoxCox(MBLi | mode ~ 1, data = blood,
     order = OR, bounds = c(0, Inf), support = c(0, 2000))
coef(as.mlt(m))
logLik(m)

## location-scale transformation model
sm <- BoxCox(MBLi ~ mode | mode, data = blood,
     order = OR, bounds = c(0, Inf), support = c(0, 2000))
coef(as.mlt(sm))
logLik(sm)

## ----STRAT-plot, fig.height = 4.6, eval = TRUE--------------------------------
## grid of delivery mode
nd <- data.frame(mode = factor(levels(blood$mode), levels = levels(blood$mode)))
q <- seq(0, max(blood$MBL), by = 1)

## plot setup
xlab <- "Measured blood loss (ml)"
ylimd <- c(0, 0.004)
  
layout(matrix(1:4, nrow = 2, byrow = FALSE))

## plot stratified model
show_legend <- FALSE
fit <- m
main <- "Stratified Transformation Model"

### distribution
h <- ylim <- c(0, 1)
type <- "l"

op <- par(mgp = c(2.5, 1, 0), mar = c(2.5, 3.5, 2, 2))
plot(seq_along(q), ylim = ylim, type = "n", ylab = "Distribution",
     main = main, xlab = xlab, lwd = 1.5, cex.main = 1.05, cex.axis = .9)
abline(h = h, lty = 2, col = "lightgrey")

#### ECDF
sapply(1:nlevels(blood$mode),
  function(i) lines(ecdf(subset(blood, mode == levels(blood$mode)[i])$MBL),
                    cex = .75, col = lcol[i]))

p <- predict(fit, newdata = nd, type = "distribution", q = q)
sapply(1:ncol(p),
  function(i) {lines(x = q, y = p[, i], type = type, col = col[i])})

if (show_legend) legend("bottomright", legend = nd$mode, col = col[seq_along(nd$mode)],
                        lty = 1, bty = "n", cex = 1)

### density
type <- "l"
main <- ""
op <- par(mgp = c(2.5, 1, 0), mar = c(4, 3.5, 0.5, 2))
plot(seq_along(q), ylim = ylimd, type = "n", ylab = "Density",
     main = main, xlab = xlab, lwd = 1.5, cex.axis = .9)
abline(h = 0, lty = 2, col = "lightgrey")

d <- predict(fit, newdata = nd, type = "density", q = q)
sapply(1:ncol(d), function(i) {lines(x = q, y = d[, i], type = type, col = col[i])})

### logLiks
text(max(q) - max(q) / 3, ylimd[2] / 10 * 9.7, labels = frmtll(logLik(fit)))

par(op)

## plot location-scale model
show_legend <- TRUE
fit <- sm
main <- toUpper(stram)

### distribution
h <- ylim <- c(0, 1)
type <- "l"

op <- par(mgp = c(2.5, 1, 0), mar = c(2.5, 3.5, 2, 2))
plot(seq_along(q), ylim = ylim, type = "n", ylab = "Distribution",
     main = main, xlab = xlab, lwd = 1.5, cex.main = 1.05, cex.axis = .9)
abline(h = h, lty = 2, col = "lightgrey")

### ECDF
sapply(1:nlevels(blood$mode),
  function(i) lines(ecdf(subset(blood, mode == levels(blood$mode)[i])$MBL),
                    cex = .75, col = lcol[i]))

p <- predict(fit, newdata = nd, type = "distribution", q = q)
sapply(1:ncol(p),
  function(i) {lines(x = q, y = p[, i], type = type, col = col[i])})

if (show_legend) legend("bottomright", legend = nd$mode, col = col[seq_along(nd$mode)],
                        lty = 1, bty = "n", cex = 1)

### density
type <- "l"
main <- ""
op <- par(mgp = c(2.5, 1, 0), mar = c(4, 3.5, 0.5, 2))
plot(seq_along(q), ylim = ylimd, type = "n", ylab = "Density",
     main = main, xlab = xlab, lwd = 1.5, cex.axis = .9)
abline(h = 0, lty = 2, col = "lightgrey")

## lines
d <- predict(fit, newdata = nd, type = "density", q = q)
sapply(1:ncol(d), function(i) {lines(x = q, y = d[, i], type = type, col = col[i])})

## logLiks
text(max(q) - max(q) / 3, ylimd[2] / 10 * 9.7, labels = frmtll(logLik(fit)))

par(op)
layout(matrix(1))


## ----STRAT_PI-----------------------------------------------------------------
PI.stram <- function(object, newdata = model.frame(sm), reference = 0) {
  stopifnot(sm$model$todistr$name == "normal")
  mu <- predict(object, type = "lp", newdata = newdata, what = "shifting")
  sigma <- sqrt(exp(-predict(object, type = "lp", newdata = newdata,
                             what = "scaling")))
  
  if (reference == 0) {
    refmu <- 0
    refsigma <- 1
  } else {
    refmu <- predict(object, type = "lp", newdata = reference,
                     what = "shifting")
    refsigma <- sqrt(exp(predict(object, type = "lp",
                                 newdata = newdata, what = "scaling")))
  }
  
  object$model$todistr$p((sigma * mu - refsigma * refmu) / 
                           sqrt(sigma^2 + refsigma^2))
}
(estPI <- PI.stram(sm, newdata = nd[nd$mode == "Cesarean section",, drop = FALSE]))

library("mvtnorm")
smm <- sm
vc <- vcov(sm)
rx <- rmvnorm(10000, coef(sm), vc)
FUN <- function(cf) {
  coef(smm) <- cf
  PI.stram(smm, newdata = nd[nd$mode == "Cesarean section",, drop = FALSE])
}
retPI <- apply(rx, MARGIN = 1, FUN = FUN)
(ciPI <- quantile(retPI, probs = c(.025, .975)))


## Section 3.1.2. Crossing hazards ##
## ----XH-setup, include = FALSE------------------------------------------------
## libraries
library("tram")
library("survival")
library("coin")
library("mpr")

## ----XH-preproc, include = FALSE----------------------------------------------
## gastric cancer data: DOI: 10.1002/1097-0142(19820501)49:9<1771::AID-CNCR2820490907>3.0.CO;2-M
data("gastric", package = "KONPsurv")

## groups
gastric$group <- factor(gastric$group, levels = c(1, 2),
                          labels = c("Control", "Intervention"))
table(gastric$group)

## Surv object
gastric$y <- with(gastric, Surv(time, status))

## grid of group
nd <- data.frame(group = factor(levels(gastric$group)))

## plot setup
ylab <- "Probability of survival"
xlab <- "Days"
q <- seq(0, 1500)
lty <- 1:3

## order of stram
OR <- 6


## ----XH-models----------------------------------------------------------------
## Weibull (mpr)
mp <- mpr(y ~ list(~ group, ~ group), data = gastric, family = "Weibull") ## always with intercept
prmp <- predict(mp, type = "survivor", newdata = nd, tvec = q)

## Cox model (tram)
sm <- Coxph(y ~ group | group, data = gastric, log_first = TRUE, order = OR)
prsm <- predict(sm, type = "survivor", newdata = nd, q = q)

## t at crossing
## t = h^{-1}[beta / (1 - sqrt(exp(gamma)))] => h^{-1}(q) = Q(F_Z(q))
sh <- predict(sm, newdata = nd[2,,drop = FALSE], type = "lp", what = c("shifting"))
scl <- sqrt(exp(predict(sm, newdata = nd[2,,drop = FALSE], type = "lp",  what = c("scaling"))))
p <- sm$model$todistr$p(sh / (1 - scl))
tx <- predict(sm, newdata = nd[2,,drop = FALSE], prob = p, type = "quantile")
sx <- predict(sm, newdata = nd[2,,drop = FALSE], q = tx, type = "survivor")

### check different parametrisations
sr1 <- Survreg(y ~ group | group, data = gastric)
sr2 <- Survreg(y ~ group | group, data = gastric, remove_intercept = FALSE)
stopifnot(isTRUE(all.equal(c(logLik(sr1)), mp$model$loglike, check.attributes = FALSE)))
stopifnot(isTRUE(all.equal(c(logLik(sr2)), mp$model$loglike, check.attributes = FALSE,
                           tol = .Machine$double.eps^{1/4})))


## ----XH-plot-model, fig.height = 3.4, fig.width = 6---------------------------
layout(matrix(1:2, nrow = 1, byrow = FALSE))

## plot Weibull model
show_legend <- FALSE
surv <- t(prmp)
main <- "Location-Scale Weibull Model"
ll <- unlist(unname(mp$model["loglike"]))
h <- ylim <- c(0, 1)
type <- "l"

op <- par(mgp = c(2.5, 1, 0), mar = c(2.5, 3.5, 2, 2), las = 1)

#### KM
plot(sf <- survfit(Surv(time, status) ~ group, data = gastric),
     ylab = "Probability of survival", main = main, xlab = xlab, 
     col = lcol, ylim = ylim, xlim = range(q),
     lwd = 1.2, cex.main = .85, cex.axis = .85)
abline(h = h, lty = 2, col = "lightgrey")

sapply(1:ncol(surv),
       function(i) {lines(x = q, y = surv[, i], type = type, col = col[i], lwd = 1.5)})

if (show_legend) legend("topright", legend = nd$group,
                        col = col[seq_along(nd$group)], lty = 1, bty = "n", cex = .8)

### logLiks
text(min(q) + 350, ylim[1] + .05, labels = frmtll(ll), cex = .85)

par(op)

## plot location-scale transformation model
show_legend <- TRUE
surv <- prsm
main <- toUpper(stram)
ll <- logLik(sm)
h <- ylim <- c(0, 1)
type <- "l"

op <- par(mgp = c(2.5, 1, 0), mar = c(2.5, 3.5, 2, 2), las = 1)

#### KM
plot(sf <- survfit(Surv(time, status) ~ group, data = gastric),
     ylab = "Probability of survival", main = main, xlab = xlab, 
     col = lcol, ylim = ylim, xlim = range(q),
     lwd = 1.2, cex.main = .85, cex.axis = .85)
abline(h = h, lty = 2, col = "lightgrey")

sapply(1:ncol(surv),
       function(i) {lines(x = q, y = surv[, i], type = type, col = col[i], lwd = 2)})

if (show_legend) legend("topright", legend = nd$group,
                        col = col[seq_along(nd$group)], lty = 1, bty = "n", cex = .8)

### logLiks
text(min(q) + 350, ylim[1] + .05, labels = frmtll(ll), cex = .85)

par(op)
layout(matrix(1))

## ----XH-test------------------------------------------------------------------
## null model & scores
m0 <- Coxph(y ~ 1, data = gastric, order = OR)

r <- resid(m0, what = "shifting")
rs <- resid(m0, what = "scaling")

### logrank test
survdiff(formula = y ~ group, data = gastric)
(plr <- pvalue(independence_test(r ~ group, data = gastric)))

### bivariate test <- logrank + scale test combined
(plrsM <- pvalue(independence_test(r + rs ~ group, data = gastric)))
(plrsQ <- pvalue(independence_test(r + rs ~ group, data = gastric, teststat = "quad")))


## Section 3.1.3 Partial proportional hazards ##
## ----PPH-setup, include = FALSE-----------------------------------------------
## libraries
library("tram")
library("cotram")
library("survival")
library("multcomp")
library("gamlss")
library("gamlss.cens")

## ----PPH-data, echo = FALSE, results = "asis", message = FALSE, warning = FALSE----
## deer-vehicle collision data: DOI: 10.1016/j.aap.2015.04.037
tmpd <- tempdir()
file <- file.path(tmpd, "analysis", "DVC.rda")
tgz <- file.path(tmpd, "DVC.tgz")
url <- "https://zenodo.org/record/17179/files"
if (!file.exists(file)) {
    op <- options(timeout = 120)
    download.file(url = paste(url, "DVC.tgz", sep = "/"), destfile = tgz)
    options(op)
    untar(tgz, exdir = tmpd)
}
load(file)

## ----PPH-preproc, echo = FALSE, results = "hide", message = FALSE-------------
## setup
loc <- Sys.setlocale("LC_ALL", "en_US.UTF-8")
rm(loc)

df <- data.frame(margin.table(obs[,,"wild"], 2))
colnames(df) <- c("day", "DVC")
df$day <- as.Date(df$day)

## baseline: Monday 2002
df$weekday <- factor(format(df$day, "%A"),
                      levels = c("Monday", "Tuesday", "Wednesday", "Thursday",
                                 "Friday", "Saturday", "Sunday"))
df$weekday[weekdays$ArbZeitFaktor == 0] <- "Sunday"
df$year <- factor(format(df$day, "%Y"))
df$time <- as.numeric(difftime(df$day, start, unit = "days"))

## harmonics
sfm <- function(timevar, freq = 365, S = 10) {
  S <- 1:S * 2
  paste("I(", rep(c("sin", "cos"), length(S)), "(",
        rep(S, rep(2, length(S))), "* pi *", timevar, "/", freq, "))",
        collapse = "+")
}

Xtime <- model.matrix(as.formula(paste("~", sfm("time"))), data = df)[,-1]
Xtime <- as.data.frame(Xtime)
colnames(Xtime) <- tvars <- paste0("tvar", 1:ncol(Xtime))

d <- cbind(df, Xtime)


## ----PPH-model, echo = FALSE, message = FALSE, warning = FALSE, results = "hide", cache = TRUE----
## formulae
lvars <- c("weekday", "year", tvars)
svars <-  tvars
mu <- paste(lvars, collapse = " + ")
sigma <- paste(svars, collapse = " + ")
fm <- as.formula(paste("DVC ~", mu)) ## location only
sfm <- as.formula(paste("DVC ~", mu, "|", sigma)) ## location-scale

## location count transformation model
m <- cotram(fm, data = d, method = "cloglog")
logLik(m)

## location-scale count transformation model
sm <- cotram(sfm, data = d, method = "cloglog")
logLik(sm)

## location-scale count Weibull model (stram)
smW <- cotram(sfm, data = d, method = "cloglog", order = 1)
logLik(smW)

## check different parametrisations
### location-scale Weibull model (gamlss with interval-censored likelihood)
gen.cens(WEI, type = "interval")
#### correct likelihood
d$DVC <- as.integer(d$DVC)
d$DVCi <- with(d, Surv(ifelse(DVC > 0, DVC - 1, -Inf), DVC, type = "interval2"))

fmi <- as.formula(paste("DVCi ~", mu))
gm <- gamlss(formula = fmi, sigma.fo = as.formula(paste("~", sigma)),
  data = d, family = WEIic, control = gamlss.control(n.cyc = 300, trace = FALSE))
### WEI2() would be the equivalent model, but does not converge
logLik(gm)
stopifnot(abs(c(logLik(smW) - logLik(gm))) < 20)

## ----PPH-grid, echo = FALSE, message = FALSE, results = "hide", warning = FALSE----
## newdata (with reference level year = "2002", and weekday = "Monday")
nd <- subset(d, year == "2002")
nd$weekday <- factor("Monday", levels = levels(nd$weekday))
nd$year <- factor("2002", levels = levels(nd$year))

q <- seq(min(df$DVC), max(df$DVC), by = .5)

## grid for levelplot
tmp <- expand.grid(DVC = q, day = nd$day)
tmp <- merge(tmp, nd[, -(which(names(nd) == "DVC"))], by = "day")
tmp$day <- as.Date(tmp$day)

## ----PPH-plot, echo = FALSE, message = FALSE, fig.height = 8.5, results='hide', out.width = ".9\\textwidth", warning = FALSE----
h <- 1:4 * 50
v <- as.numeric(as.Date(paste("2002", 1:12, "01", sep = "-")))
at <- as.Date(paste("2002", 1:12, "01", sep = "-"))
panel_DVC <- function(text, ...) {
  panel.abline(h = h, v = v, col = "lightgrey")
  panel.levelplot(...)
  panel.text(max(as.numeric(tmp$day) - 60), 190, label = text)
}

## location count transformation model
model <- m
main <- toUpper(tram)

### predict quantiles
tmp$z <- c(predict(model, type = "distribution", newdata = nd, q = q, smooth = TRUE))

### plot quantiles
loglik <- frmtll(logLik(model))
p <- levelplot(z ~ day * DVC, data = tmp, panel = panel_DVC, 
          at = 1:3 / 4, contour = TRUE, labels = list(labels = TRUE, cex = .7),
          pretty = FALSE, col.regions = rgb(.3, .3, .3, .3), colorkey = NULL,
          scales = list(x = list(at = at, format = "%b"), alternating = 1),
          main = main, ylab = "DVCs", xlab = "Day of year", text = loglik, 
          par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)
pm <- p

## location-scale count Weibull model (stram)
model <- smW
main <- "Location-Scale Weibull Model"

### predict quantiles
tmp$z <- c(predict(model, type = "distribution", newdata = nd, q = q, smooth = TRUE))

### plot quantiles
loglik <- frmtll(logLik(model))
p <- levelplot(z ~ day * DVC, data = tmp, panel = panel_DVC, 
          at = 1:3 / 4, contour = TRUE, labels = list(labels = TRUE, cex = .7),
          pretty = FALSE, col.regions = rgb(.3, .3, .3, .3), colorkey = NULL,
          scales = list(x = list(at = at, format = "%b"), alternating = 1),
          main = main, ylab = "DVCs", xlab = "Day of year", text = loglik, 
          par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)
psmW <- p

## location-scale count transformation model (stram)
model <- sm
main <- toUpper(stram)

### predict quantiles
tmp$z <- c(predict(model, type = "distribution", newdata = nd, q = q, smooth = TRUE))

### plot quantiles
loglik <- frmtll(logLik(model))
p <- levelplot(z ~ day * DVC, data = tmp, panel = panel_DVC, 
          at = 1:3 / 4, contour = TRUE, labels = list(labels = TRUE, cex = .7),
          pretty = FALSE, col.regions = rgb(.3, .3, .3, .3), colorkey = NULL,
          scales = list(x = list(at = at, format = "%b"), alternating = 1),
          main = main, ylab = "DVCs", xlab = "Day of year", text = loglik, 
          par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)
psm <- p

grid.arrange(pm, psmW, psm)


## ----PPH-seq, echo = FALSE, results = "asis", warning = FALSE-----------------
### fake contrast matrix
K <- glht(lm(DVC ~ year, data = d), mcp(year = "Sequen"))$linfct[, -1]
cf <- coef(sm)
y <- grep("year", names(cf))
cfy <- cf[y]

vc <- vcov(sm)
vcy <- vc[y, y]

### yearly multiplicative change in hazards
HR <- exp(c(1, -1)[sm$negative + 1L] * confint(glht(parm(cfy, vcy), linfct = K))$confint)
if (sm$negative) HR[,2:3] <- HR[,3:2]
HR <- formatC(HR, format = "f", digits = 2)
rownames(HR) <- sub(" - ", " -- ", rownames(HR))
tab <- data.frame(rownames(HR), cbind(HR[, 1], paste0(HR[, 2], " -- ",HR[, 3, drop = FALSE])))
colnames(tab) <- c("Year", "Hazard ratio", "95\\% CI")
print(xtable(tab, align = "llrr"), floating = FALSE, include.rownames = FALSE,
      sanitize.colnames.function = function(x){x}, 
      hline.after = c(-1, 0, nrow(tab)), booktabs = TRUE)


## Section 3.2. Location-scale transformation tree ##
## ----TRTF-setup, include = FALSE----------------------------------------------
## libraries
library("tram")
library("TH.data")
library("trtf")
library("ATR")

## ----TRTF-preproc, include = FALSE--------------------------------------------
## CHFLS data: DOI: 10.1016/j.evolhumbehav.2008.11.002.
load(file.path(path.package(package = "TH.data"), "rda", "CHFLS.rda"))

### choose necessary variables (from multcomp vignette)
org <- chfls1[, c("REGION6", "ZJ05", "ZJ06", "A35", "ZJ07", "ZJ16M", "INCRM",
                  "JK01", "JK02", "JK20", "HY04", "HY07", "A02", "AGEGAPM", 
                  "A07M", "A14", "A21", "A22M", "A23", "AX16", "INCAM", "SEXNOW", "ZW04")]

names(org) <- c("Region",
                "Rgender",               ### gender of respondent
                "Rage",                  ### age of respondent
                "RagestartA",            ### age of respondent at beginning of relationship with partner A
                "Redu",                  ### education of respondent
                "RincomeM",              ### rounded monthly income of respondent
                "RincomeComp",           ### imputed monthly income of respondent
                "Rhealth",               ### health condition respondent
                "Rheight",               ### respondent's height
                "Rhappy",                ### respondent's happiness
                "Rmartial",              ### respondent's marital status
                "RhasA",                 ### R has current A partner
                "Agender",               ### gender of partner A
                "RAagegap",              ### age gap
                "RAstartage",            ### age at marriage
                "Aheight",               ### height of partner A
                "Aedu",                  ### education of partner A
                "AincomeM",              ### rounded partner A income
                "AincomeEst",            ### estimated partner A income
                "orgasm",                ### orgasm frequency
                "AincomeComp",           ### imputed partner A income
                "Rsexnow",               ### has sex last year
                "Rhomosexual")           ### R is homosexual

### duration of partnership 
org$RAduration <- org$Rage - org$RagestartA

### code missing values
org$AincomeM[org$AincomeM < 0] <- NA
org$RincomeM[org$RincomeM < 0] <- NA
org$Aheight[org$Aheight < 0] <- NA

## outcome
olevels <- c("never", "rarely", "sometimes", "often", "always")
orgA <- subset(org, Rgender == "female" & Rhomosexual != "yes" &
                 Agender != "woman" & orgasm %in% olevels)

orgA$orgasm <- ordered(as.character(orgA$orgasm), levels = olevels)

orgA$Reducation <- factor(as.character(orgA$Redu),
        levels = rev(c("univ/grad", "j col", "up mid", "low mid", "primary", "no school")), 
        labels = rev(c("university", "junior college", "upper-middle school",
                       "lower-middle school", "primary school", "no school")), 
        ordered = TRUE)

orgA$Aeducation <- factor(as.character(orgA$Aedu),
        levels = rev(c("univ/grad", "j col", "up mid", "low mid", "primary", "no school")), 
        labels = rev(c("university", "junior college", "upper-middle school",
                       "lower-middle school", "primary school", "no school")), 
        ordered = TRUE)

orgA$Rhappy <- factor(as.character(orgA$Rhappy),
        levels = c("v unhappy", "not too", "relatively", "very"), 
        labels =  c("very unhappy", "not too unhappy", "relatively happy", "very happy"), ordered = TRUE)

orgA$Rhealth <- factor(as.character(orgA$Rhealth),
        levels = c("poor", "notgood", "fair", "good", "excellent"), ordered = TRUE)

orgA$Rregion <- factor(as.character(orgA$Region),
        levels = c("CentralW", "Northeast", "North", "InlandS", "CoastalE", "CoastalS"))

orgA$edudiff <- as.numeric(orgA$Aeducation) - as.numeric(orgA$Reducation)
orgA$wealthdiff <- orgA$RincomeComp - orgA$AincomeComp
orgA$Aincome <- orgA$AincomeComp

## ----TRTF-model---------------------------------------------------------------
## variable list taken from tram::mtram vignette
## null model
m0 <- Polr(orgasm ~ 1, data = orgA, method = "logistic")

## split wrt all baseline parameters (parm = NULL: take out baseline trafo for m0)
set.seed(20220608)
tr <- trafotree(m0, formula = orgasm ~ 1 | Aincome + Aheight + RAduration +
  Rage + edudiff + wealthdiff + Reducation + Rhealth + Rhappy + Rregion,
  data = orgA, intercept = "shift-scale", parm = NULL, maxsurrogate = 3L)
logLik(tr)


## ----TRTF-plot, fig.width = 14, fig.height = 8, out.width = "\\textwidth"-----
## this is a _very_ dirty hack, for plotting only
lev <- levels(tr$data$Rregion)
lev[1] <- paste("Other", paste(rep(" ", 100), collapse = ""))
levels(tr$data$Rregion) <- lev

plot(rotate(tr), tp_args = list(newdata = model.frame(tr)[1, -1, drop = FALSE],
  type = "density", fill = "lightgrey"), terminal_panel = trtf:::node_mlt, cex = 1.1)


## Section 3.3. Transformation additive models for location and scale ##
## not run by default (computationally too intensive)
if (run_all) {
## ----TAMLS-setup--------------------------------------------------------------
library("tram")
library("gamlss")

set.seed(2925)

## ----TAMLS-model--------------------------------------------------------------
## extract expressions for gamlss.family
e <- expression(
  d <- data.frame(y = y, m = mu, s = sigma),
  mf <- mlt(._mff, data = d, fixed = c("m" = 1, "scl_s" = 1), scale = TRUE), # theta = ._start
  if ((i/5)%%1 == 0) {print(logLik(mf))}, i <<- i + 1,
  # print(logLik(mf)),
  trm <- predict(mf, newdata = d, type = "trafo"),
  tr <- predict(mf, newdata = data.frame(y = y, s = sigma, m = 0), type = "trafo"))


## gamlss.dist for transformation models (scale_shift = FALSE)
TM <- function(mu.link = "identity", sigma.link = "identity") {
  mstats <- checklink("mu.link", "TM", substitute(mu.link), c("identity"))
  dstats <- checklink("sigma.link", "TM", substitute(sigma.link), c("identity"))
  structure(list(family = c("TM", "trafo"), 
    parameters = list(mu = TRUE, sigma = TRUE), nopar = 2, 
    type = "Continuous",
    mu.link = as.character(substitute(mu.link)), 
    sigma.link = as.character(substitute(sigma.link)),
    mu.linkfun = mstats$linkfun, 
    sigma.linkfun = dstats$linkfun,
    mu.linkinv = mstats$linkinv, 
    sigma.linkinv = dstats$linkinv,
    mu.dr = mstats$mu.eta, 
    sigma.dr = dstats$mu.eta,
    dldm = function(y, mu, sigma) {eval(e); return(trm)},
    d2ldm2 = function(y, mu, sigma) -1,
    dldd = function(y, mu, sigma) {eval(e); return(-trm * 1/2 * tr + 1/2)}, 
    d2ldd2 = function(y, mu, sigma) {eval(e); return(1/2 * (1/2 * tr * mu - tr^2))},
    d2ldmdd = function(y, mu, sigma) {eval(e); return(1/2 * tr)},
    G.dev.incr = function(y, mu, sigma) {eval(e);
      return(-2 * mf$logliki(coef(as.mlt(mf)), weights = weights(mf)))},
    rqres = expression(NA), 
    mu.initial = expression({mu <- rep(0, length(y))}),
    sigma.initial = expression({sigma <- rep(0, length(y))}),
    mu.valid = function(mu) TRUE,
    sigma.valid = function(sigma) TRUE,
    y.valid = function(y) TRUE,
    mean = function(mu, sigma) return(NA),
    variance = function(mu, sigma) return(NA)), 
  class = c("gamlss.family", "family"))
}

## get coefficients basis functions
refit <- function(model) {
  mu <- predict(model, what = "mu")
  sigma <- predict(model, what = "sigma")
  y <- dd$y
  eval(e)
  return(mf)
}

## ----TAMLS-data---------------------------------------------------------------
## head circumference data: DOI: 10.1203/00006450-200003000-00006
data("db", package = "gamlss.data")
pwr <- 0.28
db$agepwr <- with(db, age^pwr)

## ----TAMLS-model-setup--------------------------------------------------------
## support & thetas
OR <- 6
log_first <- TRUE

dd <- with(db, data.frame(y = head, m = agepwr, s = agepwr))
mf <- BoxCox(y ~ 1, data = dd, order = OR, log_first = log_first) ## thetas
._mff <- BoxCox(y ~ m | s, data = dd, model_only = TRUE, order = OR, log_first = log_first) ## support
._start <- coef(as.mlt(mf))
mlt(._mff, data = dd)

## ----TAMLS-fit----------------------------------------------------------------
## GAMLSS: from DOI: 10.18637/jss.v023.i07
mBCT <- gamlss(formula = head ~ cs(agepwr, df = 13.77), sigma.fo = ~ cs(agepwr, df = 6.05),
               nu.fo = ~ agepwr, tau.fo = ~ agepwr,
               family = BCT(), data = db)
coefAll(mBCT)
logLik(mBCT)

## TAMLS
i <- 0
mTM <- gamlss(formula = head ~ 0 + cs(agepwr, df = 13.77), sigma.fo = ~ 0 + cs(agepwr, df = 6.05),
  data = db, family = TM(), control = gamlss.control(n.cyc = 400, c.crit = 0.001))
coefAll(mTM)
logLik(mTM)

mlt_TM <- refit(mTM) ## fitted mlt model
coef(mlt_TM)
logLik(mlt_TM)

## ----TAMLS-plot, fig.height = 6.8, out.width = ".9\\textwidth", warning=FALSE, message=FALSE----
## plot setup 
pfun <- function(x, y, z, subscripts, at, text, ...) {
  panel.contourplot(x, y, z, subscripts,
    at = c(0.4, 2, 10, 25, 50, 75, 90, 98, 99.6) / 100, ...)
  panel.xyplot(x = db$age, y = db$head, pch = 20,
    col = rgb(.1, .1, .1, .1), ...)
  panel.text(15, 36, label = text)
}

## grid for head and age
nd <- expand.grid(head = sort(unique(db$head)), age = sort(unique(db$age)))
nd$agepwr <- with(nd, age^pwr)

nd$cut <- factor(nd$age > 2.5, levels = c(FALSE, TRUE),
                 labels = c("Age < 2.5 years", "Age > 2.5 years"))

## BCT GAMLSS
### quantiles
pr <- nd
pr$agepwr <- nd$agepwr <- nd$age^pwr
pr <- cbind(pr, predictAll(mBCT, newdata = nd[, c("head", "agepwr")], type = "response"))
pr$p <- with(pr, pBCT(q = y, mu = mu, sigma = sigma, nu = nu, tau = tau))

### plot
main <- "BCT GAMLSS "
loglik <- frmtll(logLik(mBCT))
p <- contourplot(p ~ age + head | cut, data = pr, panel = pfun, region = FALSE, labels = list(cex = .8),
  xlab = "Age (years)", ylab = "Head circumference (cm)", main = main, text = loglik,
  scales = list(x = list(relation = "free"), alternating = 1), layout = c(2, 1), 
  par.setting = list(layout.heights = list(bottom.padding = 0))
)
pg <- p

## TAMLS
### quantiles
pr <- nd
pr$s <- predict(mTM, what = "sigma", newdata = nd[, c("head", "agepwr")])
pr$m <- predict(mTM, what = "mu", newdata = nd[, c("head", "agepwr")])
pr$p <- c(predict(mlt_TM, newdata = with(pr, data.frame(y = head, m = m, s = s)),
  type = "distribution"))

### plot
main <- "TAMLS"
loglik <- frmtll(logLik(mTM))
## plot
p <- contourplot(p ~ age + head | cut, data = pr, panel = pfun, region = FALSE, labels = list(cex = .8),
  xlab = "Age (years)", ylab = "Head circumference (cm)", main = main, text = loglik,
  scales = list(x = list(relation = "free"), alternating = 1), layout = c(2, 1), 
  par.setting = list(layout.heights = list(bottom.padding = 0))
)
ptm <- p

grid.arrange(pg, ptm)
}

## Section 3.4. Model selection ##
## not run by default (computationally too intensive)
if (run_all) {
## ----VS-libraries-------------------------------------------------------------
## libraries
library("tram")
library("tramvs")
library("cotram")
  
set.seed(2529)
  
## ----VS-preproc, echo = FALSE-------------------------------------------------
## medical care data: https://www.jstor.org/stable/2285252
data("NMES1988", package = "AER")
nmes <- NMES1988
nmes$sex <- nmes$gender
  
mm <- model.matrix(~ health, data = nmes)[,-1]
nmes <- cbind(nmes, mm)
  
## formula
vars <- c("healthpoor", "healthexcellent", "chronic", "sex", "school",
          "insurance")
  
mu <- sigma <- paste(vars, collapse = " + ")
fm <- as.formula(paste("visits ~ " , mu)) ## location only
sfm <- as.formula(paste("visits ~ " , mu, "|", sigma)) ## location-scale
  
  
## ----VS-model, echo = FALSE, purl = TRUE, eval = FALSE------------------------
order <- 12

## transformation model maximising the count likelihood (ML)
(mML <- cotram(sfm, data = nmes, method = "cloglog",
              remove_intercept = FALSE, log_first = FALSE, 
              order = order))
  
## best subset transformation model with L0 penalty on scale term only (BSS)
(mBSS <- cotramVS(sfm, mandatory = fm, data = nmes, method = "cloglog",
                 remove_intercept = FALSE, log_first = FALSE, 
                 order = order))
mBSS <- mBSS$best_fit
                 

## ----VS-table, results = "asis"-----------------------------------------------
## ML model
vars <- c("health", vars[-grep("health", vars)]) ## model matrix
cf <- coef(mML)
cf[mML$shiftcoef] <- c(-1, 1)[mML$negative + 1L] * cf[mML$shiftcoef]
ifx <- sapply(nmes[, vars], is.factor)
f <- vars[ifx]
lev <- lapply(nmes[f], levels)
fc <- unname(do.call("c", lapply(f, function(x) {paste0(x, lev[[x]])})))
i <- c(fc, vars[!ifx])
tab <- data.frame(variable = toUpper(c(rep(f, times = unname(lengths(lev))), vars[!ifx])),
                  levels = NA, beta = cf[i], gamma = cf[paste0("scl_", i)])
tab$levels[seq_along(unlist(lev))] <- unlist(lev)
row.names(tab) <- i

tab <- tab[!is.na(tab[3]),]
colnames(tab) <- c("Variable", "Level", "$\\hat{\\beta}$","$\\hat{\\gamma}$")
tab[3:4] <- frmt4(as.matrix(tab[3:4]))
tabf <- tab

## BSS model
cf <- coef(mBSS$mod)
cf[mBSS$mod$shiftcoef] <- c(-1, 1)[mBSS$mod$negative + 1L] * cf[mBSS$mod$shiftcoef]
ifx <- sapply(nmes[, vars], is.factor)
f <- vars[ifx]
lev <- lapply(nmes[f], levels)
fc <- unname(do.call("c", lapply(f, function(x) {paste0(x, lev[[x]])})))
i <- c(fc, vars[!ifx])
tab <- data.frame(variable = toUpper(c(rep(f, times = unname(lengths(lev))), vars[!ifx])),
                  levels = NA, beta = cf[i], gamma = cf[paste0("scl_", i)])
tab$levels[seq_along(unlist(lev))] <- unlist(lev)
row.names(tab) <- i

tab <- tab[!is.na(tab[3]),]
colnames(tab) <- c("Variable", "Level", "$\\hat{\\beta}$","$\\hat{\\gamma}$")
tab[3:4] <- frmt4(as.matrix(tab[3:4]))

## cosmetics
tab <- cbind(tabf, tab[ -(1:2)]) ## combine
tab[sub("scl_", "", mBSS$I), ncol(tab)] <- "---" ## remove "unselected entries"
tab[2, 1] <- "" ## remove replicated variable name

## print table
xtab <- xtable(tab, align = "lllrrrr")
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- paste0("\\toprule &", paste0("&\\multicolumn{2}{c}{", c("ML", "BSS"), "}", collapse = ""),
                           "\\\\\\cmidrule(l){3-4} \\cmidrule(l){5-6}")
print(xtab, floating = FALSE, sanitize.colnames.function = function(x){x}, 
      sanitize.text.function = identity,
      math.style.negative = FALSE, add.to.row = addtorow, hline.after = c(0, nrow(tab)), 
      include.rownames = FALSE, booktabs = TRUE)
}


## Supplementary Material B. Re-analysis of DOI: 10.1177/0272989x8800800309 ##
## ----ROC-setup, include = FALSE-----------------------------------------------

## libraries
library("tram")

## functions ROC
ROC.stram <- function(object, newdata = model.frame(object), prob = 1:99 / 100, ...) {
  beta <- predict(object, type = "lp", newdata = newdata, what = "shifting")
  gamma <- predict(object, type = "lp", newdata = newdata, what = "scaling")
  lp <- cbind(beta, gamma)
  
  roc <- function(lp) {
    z0 <- object$model$todistr$q(1 - prob)
    beta <- lp[1]; gamma <- lp[2]
    z <- sqrt(exp(gamma)) * z0 + c(1, -1)[object$negative + 1L] * beta
    if (object$model$scale_shift) z <- sqrt(exp(gamma)) * (z0 + c(1, -1)[object$negative + 1L] * beta)
    1 - object$model$todistr$p(z)
  }
  ret <- apply(lp, 1, roc)
    
  attr(ret, "prob") <- prob
  class(ret) <- "ROCstram"
  ret
}

plot.ROCstram <- function(x, col = "black", fill = "lightgrey",
                          lty = 1, lwd = 1, add = FALSE,
                          xlab = "1 - Specificity",
                          ylab = "Sensitivity", ...) {
    prob <- attr(x, "prob")
    if (!add) {
    plot(0, 1, xlim = c(0, 1), ylim = c(0, 1), type = "n", 
         xlab = xlab, ylab = ylab, ...)
    abline(a = 0, b = 1, col = "lightgrey")
    abline(v = 0:10/10, col = "lightgrey", lty = 3)
    abline(h = 0:10/10, col = "lightgrey", lty = 3)
    }
    if (length(fill) != ncol(x)) fill <- rep(fill, length.out = ncol(x))
    if (length(col) != ncol(x)) col <- rep(col, length.out = ncol(x))
    if (length(lty) != ncol(x)) lty <- rep(lty, length.out = ncol(x))
    if (length(lwd) != ncol(x)) lwd <- rep(lwd, length.out = ncol(x))

    for (i in 1:ncol(x))
        lines(c(0, prob, 1), c(0, x[, i], 1), lty = lty[i], col = col[i], lwd = lwd[i])
}

## ----ROC-data, echo = FALSE, results = "asis", message = FALSE, warning = FALSE----
## ultrasound data: DOI: 10.1177/0272989x8800800309
tmpd <- tempdir()
url <- "https://research.fredhutch.org/content/dam/stripe/diagnostic-biomarkers-statistical-center/files"
csv <- "tostbegg2.csv"
file <- file.path(tmpd, csv)

if (!file.exists(file)) {
  op <- options(timeout = 120)
  download.file(url = paste(url, csv, sep = "/"), destfile = file)
  options(op)
}
dat <- read.csv(file)

## ----ROC-preproc--------------------------------------------------------------
## coding
dat$d <- factor(dat$d, levels = c(0, 1), labels = c("no", "yes"))
dat$y <- factor(dat$y, ordered = TRUE)
dat$type <- factor(dat$type, levels = c(0, 1), labels = c("Colon", "Breast"))

dat$x1 <- model.matrix(~ d, dat)[,-1] ## hepatitis: 0 = no, 1 = yes
dat$x2 <- model.matrix(~ type, dat)[,-1] ## primary tumor side: 0 = colon, 1 = breast
dat$x3 <- ifelse(dat$x1 == 1 & dat$x2 == 1, 1, 0) ## interaction hepatitis & breast: 0 = no, 1 = yes

## ----ROC-model----------------------------------------------------------------
## model from DOI: 10.1177/0272989x8800800309
m <- Polr(y ~ x1 + x2 + x3 | x1 + x2 + x3, data = dat, method = "probit", scale_shift = TRUE)

## ----ROC-plot, fig.height = 3.5, fig.width = 3.5, out.width = ".45\\textwidth", results = 'hold', fig.show = 'hide'----
op <- par(mgp = c(2.5, 1, 0), mar = c(3.5, 4, 0.5, 4), cex = .9)
nd <- data.frame(x1 = 1, x2 = c(0, 1), x3 = c(0, 1)) ## breast cancer no / yes
r <- ROC.stram(m, newdata = nd)
par(mgp = c(2.5, 1, 0), mar = c(4, 4, 2.5, 1), pty = "s")
plot(r, lty = 1:2, col = 1, lwd = 1.5, las = 1, cex = .5, # asp = 1,
     ylab = "True positive ratio", xlab = "False positive ratio")
legend("bottomright", legend = c("Colon", "Breast"),
       lty = 1:2, col = 1, bty = "n", lwd = 1.5)
par(op)


## Supplementary Material D. Simulation ##
## not run by default (computationally too intensive)
if (run_all) {
## ----SIM-setup----------------------------------------------------------------
## libraries
library("tram")
library("tramvs")

set.seed(2529)

## parameters
VS <- TRUE ## run tramvs

## sample size
n <- 500

## DGP
cf <- c(0, 0.5, 1.0)

Nsim <- 100

cfs <- as.matrix(expand.grid(beta = cf, gamma = cf))

negative <- FALSE ## Colr with positive shift

SIM <- function(sim) {
  
  print(sim)
  
  x <- matrix(runif(n * 2), ncol = 2)
  
  z <- rlogis(n) ## F_Z: expit
  
  if (VS) v <- matrix(runif(n * 8), ncol = 8)
  
  DGP <- function(i) {
    
    ## disjunct location and scale
    lp <- x[, 1, drop = FALSE] %*% cfs[i, "beta"] 
    sclp <- x[, 2, drop = FALSE] %*% cfs[i, "gamma"]
    
    hY <- (z + c(-1, 1)[negative + 1L] * lp) / sqrt(exp(sclp))
    
    y <- qchisq(plogis(hY), df = df <- 3) ## h^(-1) = qchisq(expit(z))
    
    d <- data.frame(y = y, x1 = x[,1], x2 = x[,2])
    
    m <- Colr(y ~ x1 | x2, data = d, bounds = c(0, Inf), scale_shift = FALSE)
    
    ## coefficients
    cfm <- c(coef(m), cfs[i, ])
    
    if (VS) {
      ## tramvs
      V <- as.data.frame(v)
      colnames(V) <- paste0("v", 1:ncol(V))
      
      d <- cbind(d, V)
      vars <- paste(colnames(d[-1]), collapse = " + ")
      mvs <- ColrVS(as.formula(paste("y ~ ",  vars, "|", vars)), data = d,
                    bounds = c(0, Inf), scale_shift = FALSE)
      
      cfvs <- mvs$best_fit$mod$coef
      cfvs <- cfvs[-grep("^Bs", names(cfvs))]
      cfvs <- c(cfvs, cfs[i, ])
      
      ## coefficients
      return(list(cfm, cfvs))
    }
    return(cfm)
  }
  ret <- lapply(1:nrow(cfs), DGP) ## fit for different coefs
  do.call("rbind", ret)
}

## ----SIM-run------------------------------------------------------------------
## parallelise for "tramvs"
if (VS) {
  library("parallel")

  RNGkind("L'Ecuyer-CMRG") ## reproducible with multiple cores
  numCores <- detectCores() - 2 ## use all available
  ret <- mclapply(1:Nsim, SIM, mc.cores = numCores)
} else {
  ret <- lapply(1:Nsim, SIM)
}

## ----SIM-preproc, include = FALSE---------------------------------------------
ret <- data.frame(do.call("rbind", ret))
d <- do.call("rbind", ret[[1]])
dvs <- do.call("rbind", ret[[2]])

## ----SIM-plot-ML, fig.height = 4, out.width = ".9\\textwidth"-----------------
d <- melt(as.data.frame(d), measure.vars = c("x1", "scl_x2"), id.vars = colnames(cfs))

d$coef <- with(d, ifelse(variable == "scl_x2", gamma, beta))
d$coef <- factor(d$coef, levels = cf, labels = paste0("coefficient = ", cf))

d$variable <- factor(d$variable, levels = c("scl_x2", "x1"))

ylim <- c(-1.5, 2.5)

## plot
lattice.options(layout.heights = list(bottom.padding = list(x = 0), top.padding = list(x = 0)))
p1 <- bwplot(value ~ factor(gamma) | coef, data = d, subset = d$variable == "x1",
  ylab = bquote(widehat(beta)), xlab = expression(gamma), ylim = ylim, layout = c(3, 1),
  scales = list(x = list(relation = "free"), alternating = 1),
  strip = strip.custom(factor.levels = c(expression(beta==0), expression(beta==0.5), expression(beta==1))),
  panel = function(x, y, groups, subscripts, ...) {
    panel.abline(h = cf[unique(as.numeric(d$coef[subscripts]))], col = "gray70", lty = 2)
    panel.bwplot(x = x, y = y, ...)},
   par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)

p2 <- bwplot(value ~ factor(beta) | coef, data = d, subset = d$variable == "scl_x2",
  ylab = bquote(widehat(gamma)), xlab = expression(beta),  ylim = ylim, layout = c(3, 1),
  scales = list(x = list(relation = "free"), alternating = 1),
  strip = strip.custom(factor.levels = c(expression(gamma==0), expression(gamma==0.5), expression(gamma==1))),
  panel = function(x, y, groups, subscripts, ...) {
    panel.abline(h = cf[unique(as.numeric(d$coef[subscripts]))], col = "gray70", lty = 2)
    panel.bwplot(x = x, y = y, ...)},
   par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)

## combine all plots 
grid.arrange(p1, p2, ncol = 1)


## ----SIM-plot-VS, fig.height = 8, out.width = ".9\\textwidth"-----------------
## informative covariates
d <- melt(as.data.frame(dvs), measure.vars = c("x1", "scl_x2"), id.vars = colnames(cfs))
d$coef <- with(d, ifelse(variable == "scl_x2", gamma, beta))
d$coef <- factor(d$coef, levels = cf, labels = paste0("coefficient = ", cf))

d$variable <- factor(d$variable, levels = c("scl_x2", "x1"))

ylim <- c(-1.5, 2.5)

## plot
lattice.options(layout.heights = list(bottom.padding = list(x = 0), top.padding = list(x = 0)))
p1 <- bwplot(value ~ factor(gamma) | coef, data = d, subset = d$variable == "x1",
  ylab = bquote(widehat(beta)), xlab = expression(gamma), ylim = ylim, layout = c(3, 1),
  scales = list(x = list(relation = "free"), alternating = 1),
  strip = strip.custom(factor.levels = c(expression(beta==0), expression(beta==0.5), expression(beta==1))),
  panel = function(x, y, groups, subscripts, ...) {
    panel.abline(h = cf[unique(as.numeric(d$coef[subscripts]))], col = "gray70", lty = 2)
    panel.bwplot(x = x, y = y, ...)},
   par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)

p2 <- bwplot(value ~ factor(beta) | coef, data = d, subset = d$variable == "scl_x2",
  ylab = bquote(widehat(gamma)), xlab = expression(beta),  ylim = ylim, layout = c(3, 1),
  scales = list(x = list(relation = "free"), alternating = 1),
  strip = strip.custom(factor.levels = c(expression(gamma==0), expression(gamma==0.5), expression(gamma==1))),
  panel = function(x, y, groups, subscripts, ...) {
    panel.abline(h = cf[unique(as.numeric(d$coef[subscripts]))], col = "gray70", lty = 2)
    panel.bwplot(x = x, y = y, ...)},
   par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)

## uninformative covariates
v <- paste0("v", 1:8)
vscl <- paste0("scl_v", 1:8)
dv <- melt(as.data.frame(dvs), measure.vars = c(v, vscl), id.vars = colnames(cfs))

dv$variable <- factor(dv$variable, levels = c(v, vscl))

dv$gamma <- factor(dv$gamma, levels = gamma <- cf, labels = paste0("gamma = ", gamma))
dv$beta <- factor(dv$beta, levels = beta <- rev(cf), labels = paste0("beta = ", beta))

i <- rep(3:10, 2)
l <- rep(c("location", "scale"), each = 8)
label <- sapply(seq_along(i), function(t) as.expression(bquote(.(l[t])~x[.(i[t])])))

## plot
trellis.par.set(list(axis.components = list(top = list(tck = 0))))
p <- xyplot(value ~ variable | gamma * beta, data = dv, 
  ylab = "Estimate", xlab = "Covariate", ylim = ylim,
  scales = list(x = list(rot = 45, cex = .6, 
    label = label), alternating = 1),
  main = "Uninformative covariates", ## better formatting / title?
  axis.components = list(top = list(tck = 0)),
  panel = function(x, y, groups, subscripts, ...) {
    panel.abline(h = 0, col = "gray70", lty = 2)
    panel.abline(v = 8.5, col = 1, lty = 2)
    panel.xyplot(x = x, y = y, pch = 19, col = rgb(.1, .1, .1, .1), ...)},
   par.settings = list(layout.heights = list(main.key.padding = 0, bottom.padding = 0))
)

pv <- useOuterStrips(p, 
  strip = strip.custom(which.given = 1,
    factor.levels = sapply(gamma, function(g) as.expression(bquote(gamma~"="~.(g))))),
  strip.left = strip.custom(which.given = 1,
    factor.levels = sapply(beta, function(b) as.expression(bquote(beta~"="~.(b))))))

## combine all plots 
grid.arrange(p1, p2, pv, heights = c(2, 2, 4))
}


## Session info ##
warnings()
sessionInfo()
