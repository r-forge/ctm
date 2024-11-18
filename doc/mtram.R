## ----mtram-pkgs, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
set.seed(290875)

pkgs <- c("colorspace", "survival", "lme4", "tram", "gridExtra",
          "lattice", "latticeExtra", "mvtnorm", "ordinalCont", "tramME")
pkgs <- sapply(pkgs, require, character.only = TRUE)


## ----mtram-citation, echo = FALSE---------------------------------------------
year <- substr(packageDescription("tram")$Date, 1, 4)
version <- packageDescription("tram")$Version


## ----fail, results = "asis", echo = FALSE-------------------------------------
if (any(!pkgs))
{
    cat(paste("Package(s)", paste(names(pkgs)[!pkgs], collapse = ", "), 
        "not available, stop processing.",
        "\\end{document}\n"))
    knitr::knit_exit()
}
if (!interactive() && .Platform$OS.type != "unix")
{
    cat("Vignette only compiled under Unix alikes.")
    knitr::knit_exit()
}


## ----mtram-setup, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
trellis.par.set(list(plot.symbol = list(col=1,pch=20, cex=0.7),
                     box.rectangle = list(col=1),
                     box.umbrella = list(lty=1, col=1),
                     strip.background = list(col = "white")))
ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)

knitr::opts_chunk$set(echo = TRUE, results = 'markup', error = FALSE,
                      warning = FALSE, message = FALSE,
                      tidy = FALSE, cache = FALSE, size = "small",
                      fig.width = 6, fig.height = 4, fig.align = "center",
                      out.width = NULL, ###'.6\\linewidth', 
                      out.height = NULL,
                      fig.scap = NA)
knitr::render_sweave()  # use Sweave environments
knitr::set_header(highlight = '')  # do not \usepackage{Sweave}
## R settings
options(prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)  # JSS style
options(width = 75)

### ecdf plots
myprepanel <- function (x, y, f.value = NULL, ...) 
{
    ans <- prepanel.default.qqmath(x, f.value = f.value, distribution = qunif)
    with(ans, list(xlim = ylim, ylim = c(0, 1), dx = dy, dy = dx))
}


mypanel <- function (x, y, f.value = NULL, type = "s", groups = NULL, qtype = 7, 
    ref = TRUE, ...) 
{
    if (ref) {
        reference.line <- trellis.par.get("reference.line")
        do.call(panel.abline, c(list(h = c(0, 1)), reference.line))
    }
    x <- as.numeric(x)
    distribution <- qunif
    nobs <- sum(!is.na(x))
    if (!is.null(groups)) {
        panel.superpose(x, y = NULL, f.value = f.value, type = type, 
            distribution = distribution, qtype = qtype, groups = groups, 
            panel.groups = panel.ecdfplot, ...)
    }
    else if (nobs) {
        if (is.null(f.value)) {
            panel.xyplot(x = sort(x), y = cumsum(y[order(x)]) / sum(y),
                type = type, ...)
        }
        else {
            p <- if (is.numeric(f.value)) 
                f.value
            else f.value(nobs)
            panel.xyplot(x = quantile(x, p, names = FALSE, type = qtype, 
                na.rm = TRUE), y = distribution(p), type = type, 
                ...)
        }
    }
}
col <- diverge_hcl(2, h = c(246, 40), c = 120, l = c(65, 90), alpha = .75)


## ----mtram-vignette, eval = FALSE---------------------------------------------
# install.packages("tram")
# demo("mtram", package = "tram")


## ----mtram-sleep-plot, echo = FALSE-------------------------------------------
library("lme4")
xyplot(Reaction ~ Days | Subject, data = sleepstudy, 
       xlab = "Days of sleep deprivation", ylab = "Average reaction time (in ms)")


## ----mtram-sleep_lmer, cache = FALSE------------------------------------------
sleep_lmer <- lmer(Reaction ~ Days + (Days | Subject), 
                   data = sleepstudy, REML = FALSE)


## ----mtram-tram---------------------------------------------------------------
library("tram")

## ----mtram-sleep_mtram, cache = FALSE-----------------------------------------
sleep_LM <- Lm(Reaction ~ Days, data = sleepstudy)
sleep_LMmer <- mtram(sleep_LM, ~ (Days | Subject), data = sleepstudy)


## ----mtram-sleep_cmp----------------------------------------------------------
logLik(sleep_lmer)
logLik(sleep_LMmer)


## ----mtram-sleep_sd-----------------------------------------------------------
(sdinv <- 1 / summary(sleep_lmer)$sigma)
coef(sleep_LMmer)["Reaction"]


## ----mtram-sleep_beta---------------------------------------------------------
fixef(sleep_lmer) * c(-1, 1) * sdinv
coef(sleep_LMmer)[c("(Intercept)", "Days")]


## ----mtram-sleep_varparm------------------------------------------------------
sleep_lmer@theta
coef(sleep_LMmer)[-(1:3)]


## ----mtram-sleep_varcorr------------------------------------------------------
sleep_LMmer$G * (1 / sdinv)^2
cov2cor(sleep_LMmer$G * (1 / sdinv)^2)
unclass(VarCorr(sleep_lmer))$Subject


## ----mtram-sleep-Surv---------------------------------------------------------
library("survival")
sleepstudy$Reaction_I <- with(sleepstudy, Surv(Reaction - 20, Reaction + 20, 
                                               type = "interval2"))
sleepstudy$Reaction_I[1:5]


## ----mtram-sleep-interval, cache = FALSE--------------------------------------
sleep_LM_I <- Lm(Reaction_I ~ Days, data = sleepstudy)
sleep_LMmer_I <- mtram(sleep_LM_I, ~ (Days | Subject), data = sleepstudy)


## ----mtram-sleep-interval-results---------------------------------------------
logLik(sleep_LMmer_I)
coef(sleep_LMmer_I)
coef(sleep_LMmer)


## ----mtram-sleep_BoxCox, cache = FALSE----------------------------------------
sleep_BC <- BoxCox(Reaction ~ Days, data = sleepstudy)
sleep_BCmer <- mtram(sleep_BC, ~ (Days | Subject), data = sleepstudy)
logLik(sleep_BCmer)


## ----mtram-sleep_BoxCoxPlot, echo = FALSE, fig.height = 5---------------------
tmp <- as.mlt(sleep_BC)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_BCmer)[names(cf)]
plot(tmp, newdata = data.frame(Days = 0), type = "trafo", col = "black",
     xlab = "Average reaction time (in ms)", ylab = expression(h(y)))
rug(sleepstudy$Reaction, col = rgb(.1, .1, .1, .1))


## ----mtram-sleep_marginal, fig.height = 5, fig.width = 7, echo = FALSE--------
days <- 0:9
q <- seq(from = min(sleepstudy$Reaction), to = max(sleepstudy$Reaction), 
         length.out = 100)
meth <- c("Normal linear mixed-effects model", "Non-normal linear transformation model")
ex <- expand.grid(Reaction = q, Days = days,
                  Method = factor(meth, levels = meth, labels = meth))
U <- cbind(1, days)
### Linear
tmp <- as.mlt(sleep_LM)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_LMmer)[names(cf)]
SLM <- tcrossprod(U %*% as(sleep_BCmer$G, "matrix"), U) + diag(nrow(U))
sd <- sqrt(diag(SLM))
h <- predict(tmp, newdata = data.frame(Days = days), q = q, type = "trafo")
prob_LM <- pnorm(t(t(h) / sd ))
### BoxCox
tmp <- as.mlt(sleep_BC)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_BCmer)[names(cf)]
SBC <- tcrossprod(U %*% as(sleep_BCmer$G, "matrix"), U) + diag(nrow(U))
sd <- sqrt(diag(SBC))
h <- predict(tmp, newdata = data.frame(Days = days), q = q, type = "trafo")
prob_BC <- pnorm(t(t(h) / sd ))
ex$prob <- c(prob_LM, prob_BC)
plotfun <- function(prob, data, ...) {
    fm <- as.formula(paste(prob, "~ Reaction | Days"))
    xyplot(fm, data = data, type = "l", 
        panel = function(x, y, subscripts, ...) {
            tmp <- subset(sleepstudy, Days == unique(nd[subscripts, "Days"]))
            mypanel(tmp$Reaction, rep(1, nrow(tmp)), lwd = 3, col = grey)
            panel.xyplot(x, y, subscripts = subscripts, ...)
    }, col = col,  xlab = "Average reaction time (in ms)", 
    ylab = "Marginal distribution function", lwd = 2, groups = Method, ...)
}
grey <- rgb(.75, .75, .75)
nd <- ex
plotfun("prob", ex, layout = c(5, 2), par.settings = simpleTheme(col=c(grey, col), lwd = 3),
  auto.key = list(text = c("Empirical cumulative distribution function", levels(nd$Method)), 
                  points = FALSE, lines = TRUE, space = "top"))


## ----mtram-sleep_corr---------------------------------------------------------
cov2cor(sleep_BCmer$G)


## ----mtram-sleep_vcov---------------------------------------------------------
library("mvtnorm")
VC <- vcov(sleep_BCmer)
idx <- (nrow(VC) - 2):nrow(VC)
Rcoef <- rmvnorm(1000, mean = coef(sleep_BCmer), sigma = VC)[,idx]
ret <- apply(Rcoef, 1, function(gamma) {
    L <- matrix(c(gamma[1:2], 0, gamma[3]), nrow = 2)
    V <- tcrossprod(L)
    c(diag(V), cov2cor(V)[1,2])
})


## ----mtram-sleep_ci-----------------------------------------------------------
### variance random intercept
quantile(ret[1,], c(.025, .5, .975))
### variance random slope
quantile(ret[2,], c(.025, .5, .975))
### correlation random intercept / random slope
quantile(ret[3,], c(.025, .5, .975))


## ----mtram-sleep_Colr, cache = FALSE------------------------------------------
sleep_C <- Colr(Reaction ~ Days, data = sleepstudy)
sleep_Cmer <- mtram(sleep_C, ~ (Days | Subject), data = sleepstudy)
logLik(sleep_Cmer)


## ----mtram-sleep_marginal-2, fig.height = 5, fig.width = 7, echo = FALSE------
days <- 0:9
q <- seq(from = min(sleepstudy$Reaction), to = max(sleepstudy$Reaction), 
         length.out = 100)
meth <- c("Normal linear mixed-effects model", "Probit transformation model", 
          "Marginal logit transformation model")
ex <- expand.grid(Reaction = q, Days = days,
                  Method = factor(meth, levels = meth, labels = meth))
U <- cbind(1, days)
### Linear
tmp <- as.mlt(sleep_LM)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_LMmer)[names(cf)]
SLM <- tcrossprod(U %*% as(sleep_BCmer$G, "matrix"), U) + diag(nrow(U))
sd <- sqrt(diag(SLM))
h <- predict(tmp, newdata = data.frame(Days = days), q = q, type = "trafo")
prob_LM <- pnorm(t(t(h) / sd ))
### BoxCox
tmp <- as.mlt(sleep_BC)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_BCmer)[names(cf)]
SBC <- tcrossprod(U %*% as(sleep_BCmer$G, "matrix"), U) + diag(nrow(U))
sd <- sqrt(diag(SBC))
h <- predict(tmp, newdata = data.frame(Days = days), q = q, type = "trafo")
prob_BC <- pnorm(t(t(h) / sd ))
### Colr
tmp <- as.mlt(sleep_C)
cf <- coef(tmp)
coef(tmp) <- coef(sleep_Cmer)[names(cf)]
SBC <- tcrossprod(U %*% as(sleep_Cmer$G, "matrix"), U) + diag(nrow(U))
sd <- sqrt(diag(SBC))
h <- predict(tmp, newdata = data.frame(Days = days), q = q, type = "trafo")
prob_C <- plogis(t(t(h) / sd ))
ex$prob <- c(prob_LM, prob_BC, prob_C)
plotfun <- function(prob, data, ...) {
    fm <- as.formula(paste(prob, "~ Reaction | Days"))
    xyplot(fm, data = data, type = "l", 
        panel = function(x, y, subscripts, ...) {
            tmp <- subset(sleepstudy, Days == unique(nd[subscripts, "Days"]))
            mypanel(tmp$Reaction, rep(1, nrow(tmp)), lwd = 3, col = grey)
            panel.xyplot(x, y, subscripts = subscripts, ...)
    }, col = c(col, col[2]),  xlab = "Average reaction time (in ms)", 
    ylab = "Marginal distribution function", lwd = 2, groups = Method, lty =
    c(1, 1, 3), ...)
}
grey <- rgb(.75, .75, .75)
nd <- ex
plotfun("prob", ex, layout = c(5, 2), par.settings = simpleTheme(col=c(grey, col, col[2]), lwd =
3, lty = c(1, 1, 1, 3)),
  auto.key = list(text = c("Empirical cumulative distribution function", levels(nd$Method)), 
                  points = FALSE, lines = TRUE, space = "top"))


## ----mtram-toenail-plot, echo = FALSE, cache = FALSE--------------------------
data("toenail", package = "HSAUR3")
layout(matrix(1:2, ncol = 2))
trt <- levels(toenail$treatment)
cdplot(outcome ~ time, data = subset(toenail, treatment == trt[1]),
       main = trt[1], xlab = "Time", ylab = "Toe nail infection")
cdplot(outcome ~ time, data = subset(toenail, treatment == trt[2]),
       main = trt[2], xlab = "Time", ylab = "")


## ----mtram-toenail_glmer_RI, cache = FALSE------------------------------------
### Laplace
toenail_glmer_RI_1 <- 
    glmer(outcome ~ treatment * time + (1 | patientID),
          data = toenail, family = binomial(link = "probit"), 
          nAGQ = 1)
summary(toenail_glmer_RI_1)
toenail_glmer_RI_1@theta

### Adaptive Gaussian Quadrature
toenail_glmer_RI_2 <- 
    glmer(outcome ~ treatment * time + (1 | patientID),
          data = toenail, family = binomial(link = "probit"), 
          nAGQ = 20)
summary(toenail_glmer_RI_2)
toenail_glmer_RI_2@theta


## ----mtram-toenail_glmmTMB_RI, cache = FALSE----------------------------------
library("glmmTMB")
toenail_glmmTMB_RI_3 <- 
    glmmTMB(outcome ~ treatment * time + (1 | patientID),
         data = toenail, family = binomial(link = "probit"))
summary(toenail_glmmTMB_RI_3)


## ----mtram-toenail_mtram_RI, cache = FALSE------------------------------------
m <- ctm(as.basis(~ outcome, data = toenail), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Normal", negative = TRUE)
toenail_probit <- mlt(m, data = toenail, 
                      fixed = c("outcomemoderate or severe" = 0))
toenail_mtram_RI <- 
    mtram(toenail_probit, ~ (1 | patientID), data = toenail)
coef(toenail_mtram_RI)


## ----mtram-toenail-hessian----------------------------------------------------
vcov(toenail_glmer_RI_2)
vcov(toenail_mtram_RI)[1:4, 1:4]


## ----mtram-toenail-coef-------------------------------------------------------
cf <- coef(toenail_mtram_RI)
cf[2:4] / sqrt(1 + cf["gamma1"]^2)


## ----mtram-toenail-gee-probit-------------------------------------------------
library("geepack")
gin <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, corstr = "independence", 
              family = binomial(link = "probit"))
gex <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, cor = "exchangeable", 
              family = binomial(link = "probit"))
gun <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, cor = "unstructured", 
              family = binomial(link = "probit"))


## ----mtram-toenail-gee-probit-coef--------------------------------------------
cbind(mtram = cf[2:4] / sqrt(1 + cf["gamma1"]^2),
      indep = coef(gin)[-1],
      excha = coef(gex)[-1],
      unstr = coef(gun)[-1])


## ----mtram-toenail-gee-logit--------------------------------------------------
gin <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, corstr = "independence", 
              family = binomial())
gex <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, cor = "exchangeable", 
              family = binomial())
gun <- geeglm(I((0:1)[outcome]) ~ treatment * time, 
              id = patientID, data = toenail, cor = "unstructured", 
              family = binomial())


## ----mtram-toenail-gee-logit-coef---------------------------------------------
coef(gin)
coef(gex)
coef(gun)


## ----mtram-toenail_logit, cache = FALSE---------------------------------------
m <- ctm(as.basis(~ outcome, data = toenail), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Logistic", negative = TRUE)
toenail_logit <- mlt(m, data = toenail, 
                     fixed = c("outcomemoderate or severe" = 0))
toenail_mtram_logit <- mtram(toenail_logit, ~ (1 | patientID), 
                             data = toenail)


## ----mtram-toenail_logit-coef-------------------------------------------------
cf <- coef(toenail_mtram_logit)
cf[2:4] / sqrt(1 + cf["gamma1"]^2)


## ----mtram-toenail-trt--------------------------------------------------------
S <- rmvnorm(10000, mean = coef(toenail_mtram_logit), 
             sigma = vcov(toenail_mtram_logit))
(ci <- quantile(S[,"treatmentterbinafine:time"] / sqrt(1 + S[, "gamma1"]^2), 
                prob = c(.025, .975)))


## ----mtram-toenail-gee-logit-mcoef--------------------------------------------
cbind(mtram = cf[2:4] / sqrt(1 + cf["gamma1"]^2),
      indep = coef(gin)[-1],
      excha = coef(gex)[-1],
      unstr = coef(gun)[-1])


## ----mtram-GEE-CI-------------------------------------------------------------
exp(coef(gun)["treatmentterbinafine:time"] +
    c(-1, 1) * qnorm(.975) * sqrt(diag(vcov(gun)))["treatmentterbinafine:time"])


## ----mtram-toenail_marginal_logit_s-------------------------------------------
tmp <- toenail_logit
cf <- coef(tmp)
cf <- cf[names(cf) != "outcomemoderate or severe"]
sdrf <- rev(coef(toenail_mtram_logit))[1]
cf <- coef(toenail_mtram_logit)[names(cf)] / sqrt(sdrf^2 + 1)
cf <- c(cf[1], "outcomemoderate or severe" = 0, cf[-1])
coef(tmp) <- cf
time <- 0:180/10
treatment <- sort(unique(toenail$treatment))
nd <- expand.grid(time = time, treatment = treatment)
nd$prob_logit <- predict(tmp, newdata = nd, type = "distribution")[1,]
nd$odds <- exp(predict(tmp, newdata = nd, type = "trafo")[1,])


## ----mtram-toenail_OR_2, dev = "png", cache = FALSE, echo = FALSE, dpi = 300----
X <- model.matrix(~ treatment * time, data = nd)
rbeta <- rmvnorm(10000, mean = coef(toenail_mtram_logit), 
                 sigma = vcov(toenail_mtram_logit))
s <- rbeta[,ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
odds <- exp(-X %*% t(rbeta))
OR <- odds[1:length(time),] / odds[-(1:length(time)),]
plot(time, rep(0, length(time)), ylim = range(OR), type = "n", 
     xlab = "Time", ylab = "Odds ratio")
colgrey <- rgb(.1, .1, .1, .01)
out <- apply(OR, 2, function(x) lines(time, x, col = colgrey))
ORest <- nd$odds[1:length(time)] / nd$odds[-(1:length(time))]
lines(time, ORest, col = col[1], lwd = 2)


## ----mtram-toenail_marginal_logit---------------------------------------------
tmp <- toenail_logit
cf <- coef(tmp)
cf <- cf[names(cf) != "outcomemoderate or severe"]
sdrf <- rev(coef(toenail_mtram_logit))[1]
cf <- coef(toenail_mtram_logit)[names(cf)] 
cf <- c(cf[1], "outcomemoderate or severe" = 0, cf[-1])
coef(tmp) <- cf
pr <- predict(tmp, newdata = nd, type = "distribution")[1,]
nd$prob_logit <- pnorm(qnorm(pr) / sdrf)


## ----mtram-toenail_marginal_probit--------------------------------------------
tmp <- toenail_probit
cf <- coef(tmp)
cf <- cf[names(cf) != "outcomemoderate or severe"]
sdrf <- rev(coef(toenail_mtram_RI))[1]
cf <- coef(toenail_mtram_RI)[names(cf)] / sqrt(sdrf^2 + 1)
cf <- c(cf[1], "outcomemoderate or severe" = 0, cf[-1])
coef(tmp) <- cf
nd$prob_probit <- predict(tmp, newdata = nd, type = "distribution")[1,]


## ----mtram-toenail_probplot, echo = FALSE-------------------------------------
nd2 <- nd[rep(1:nrow(nd), 2),]
nd2$prob <- c(nd$prob_probit, nd$prob_logit)
lev <- c("Probit", "Logit")
nd2$model <- rep(factor(lev, labels = lev, levels = lev), each = nrow(nd))

xyplot(prob ~ time | model, data = nd2, group = treatment, ylim = c(0, 1), 
       xlab = "Time", 
       par.settings = simpleTheme(col = col),
       auto.key = list(text = levels(nd2$treatment), 
                       points = FALSE, lines = TRUE, space = "top"), 
       col = col, type = "l", ylab = "Probability (none or mild)")


## ----mtram-toenail-subset-----------------------------------------------------
(rlev <- levels(toenail$patientID)[xtabs(~ patientID, 
                                        data = toenail) == 1])
toenail_gr1 <- subset(toenail, !patientID %in% rlev)
toenail_gr1$patientID <- toenail_gr1$patientID[, drop = TRUE]


## ----mtram-toenail_glmer_RS, cache = FALSE------------------------------------
toenail_glmer_RS <- 
    glmer(outcome ~ treatment * time + (1 + time | patientID),
          data = toenail_gr1, family = binomial(link = "probit"))
summary(toenail_glmer_RS)
toenail_glmer_RS@theta


## ----mtram-toenail_glmmTMB_RS, cache = FALSE----------------------------------
toenail_glmmTMB_RS_1 <- 
    glmmTMB(outcome ~ treatment * time + (1 + time | patientID),
         data = toenail_gr1, family = binomial(link = "probit"))
summary(toenail_glmmTMB_RS_1)


## ----mtram-toenail_mtram_RS, cache = FALSE------------------------------------
m <- ctm(as.basis(~ outcome, data = toenail_gr1), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Normal", negative = TRUE)
toenail_probit <- mlt(m, data = toenail_gr1, 
                      fixed = c("outcomemoderate or severe" = 0))
toenail_mtram_RS <- 
    mtram(toenail_probit, ~ (1 + time | patientID), 
          data = toenail_gr1)
logLik(toenail_mtram_RS)
coef(toenail_mtram_RS)


## ----toenail-comparisons, cache = FALSE, echo = FALSE, results = "hide"-------
t1 <- system.time(toenail_glmer_RI_1 <- 
    glmer(outcome ~ treatment * time + (1 | patientID),
          data = toenail, family = binomial(link = "probit"), 
          nAGQ = 1))

t2 <- system.time(toenail_glmer_RI_2 <- 
    glmer(outcome ~ treatment * time + (1 | patientID),
          data = toenail, family = binomial(link = "probit"), 
          nAGQ = 20))

t3 <- system.time(toenail_glmmTMB_RI_3 <- 
    glmmTMB(outcome ~ treatment * time + (1 | patientID),
         data = toenail, family = binomial(link = "probit")))

m <- ctm(as.basis(~ outcome, data = toenail), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Normal", negative = TRUE)
toenail_probit <- mlt(m, data = toenail, 
                      fixed = c("outcomemoderate or severe" = 0))
t4 <- system.time(toenail_mtram_RI <- 
    mtram(toenail_probit, ~ (1 | patientID), data = toenail))

t5 <- system.time(toenail_glmer_RS <- 
    glmer(outcome ~ treatment * time + (1 + time | patientID),
          data = toenail_gr1, family = binomial(link = "probit")))

t6 <- system.time(toenail_glmmTMB_RS_1 <- 
    glmmTMB(outcome ~ treatment * time + (1 + time | patientID),
         data = toenail_gr1, family = binomial(link = "probit")))

m <- ctm(as.basis(~ outcome, data = toenail_gr1), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Normal", negative = TRUE)
toenail_probit <- mlt(m, data = toenail_gr1, 
                      fixed = c("outcomemoderate or severe" = 0))
t7 <- system.time(toenail_mtram_RS <- 
    mtram(toenail_probit, ~ (1 + time | patientID), 
           data = toenail_gr1))

## ----output, echo = FALSE------------------------------------------------
tn_RI_glmer_L <- c(fixef(toenail_glmer_RI_1), toenail_glmer_RI_1@theta, 0, 0)
tn_RI_glmer_A <- c(fixef(toenail_glmer_RI_2), toenail_glmer_RI_2@theta, 0, 0)
tn_RI_glmmTMB <- c(fixef(toenail_glmmTMB_RI_3)$cond, sqrt(VarCorr(toenail_glmmTMB_RI_3)$cond$patientID), 0, 0)
tn_RI_mlt <- c(coef(toenail_mtram_RI), 0, 0)
tn_RS_glmer <- c(fixef(toenail_glmer_RS), toenail_glmer_RS@theta)
tn_RS_glmmTMB <- c(fixef(toenail_glmer_RS), chol(VarCorr(toenail_glmmTMB_RS_1)$cond$patientID)[c(1,3, 4)])
tn_RS_mlt <- coef(toenail_mtram_RS)
tn <- cbind(tn_RI_glmer_L, tn_RI_glmer_A , tn_RI_glmmTMB, tn_RI_mlt ,
            tn_RS_glmer, tn_RS_glmmTMB, tn_RS_mlt)

logLik(toenail_glmer_RI_1)
logLik(toenail_glmer_RI_2)
logLik(toenail_glmmTMB_RI_3)
logLik(toenail_mtram_RI)

logLik(toenail_glmer_RS)
logLik(toenail_glmmTMB_RS_1)
logLik(toenail_mtram_RS)

ll <- c(
### logLik of transformation model for glmer (Laplace) parameters
logLik(toenail_mtram_RI, tn_RI_glmer_L[1:5] * c(-1, 1, 1, 1, 1)),
### logLik of transformation model for glmer (AGQ) parameters
logLik(toenail_mtram_RI, tn_RI_glmer_A[1:5] * c(-1, 1, 1, 1, 1)),
### logLik of transformation model for glmmTMB (Laplace) parameters
logLik(toenail_mtram_RI, tn_RI_glmmTMB[1:5] * c(-1, 1, 1, 1, 1)),
### logLik of transformation model
logLik(toenail_mtram_RI),
### logLik of transformation model for glmer (Laplace) parameters
logLik(toenail_mtram_RS, tn_RS_glmer * c(-1, rep(1, 6))),
### logLik of transformation model for glmmTMB (Laplace) parameters
logLik(toenail_mtram_RS, tn_RS_glmmTMB * c(-1, rep(1, 6))),
### logLik of transformation model
logLik(toenail_mtram_RS))

tm <- c(t1["user.self"],
        t2["user.self"],
        t3["user.self"],
        t4["user.self"],
        t5["user.self"],
        t6["user.self"],
        t7["user.self"])
tm <- formatC(tm, format = "f", digits = 2, width = 5)

tn <- formatC(tn, format = "f", digits = 2, width = 5)
ll <- formatC(ll, format = "f", digits = 2, width = 6)
tn <- cbind(c("$\\alpha$", "$\\eshiftparm_1$", "$\\eshiftparm_2$", "$\\eshiftparm_3$", "$\\gamma_1$", "$\\gamma_2$", "$\\gamma_3$"), tn)
ret <- c("
\\begin{tabular}{lrrrr|rrr} \\\\ \\hline
& \\multicolumn{4}{c|}{RI} & \\multicolumn{3}{c}{RI + RS} \\\\
& \\texttt{glmer} & \\texttt{glmer} & \\texttt{glmmTMB} &  & \\texttt{glmer} & \\texttt{glmmTMB} & \\\\
& L               & AGQ             & L & (7) & L & L & (7) \\\\ \\hline")
ret <- c(ret, apply(tn, 1, function(x) c(paste(x, collapse = " & "), "\\\\")))
ret <- c(ret, "\\hline")
ret <- c(ret, 
         paste("LogLik &", paste(ll, collapse = "&"), "\\\\ "), 
         paste("Time (sec)   &", paste(tm, collapse = "&"), "\\\\ \\hline"), 
         "\\end{tabular}")

## ----table, echo = FALSE, results = "asis"------------------------------------
cat(ret, sep = "\n")


## ----mtram-neck_plot, echo = FALSE, fig.height = 4, fig.width = 7-------------
data("neck_pain", package = "ordinalCont")
pain_df <- neck_pain
idt <- xtabs(~ id, data = pain_df)
miss <- names(idt)[idt < 3]
pain_df <- subset(pain_df, !id %in% miss)
pain_df$id <- factor(pain_df$id)
levels(pain_df$laser) <- c("Placebo", "Active")
levels(pain_df$time) <- c("Baseline", "7 weeks", "12 weeks")
pain <- rbind(subset(pain_df, laser == levels(pain_df$laser)[1]),
              subset(pain_df, laser == levels(pain_df$laser)[2]))
p1 <- xyplot(vas ~ time | laser, data = pain, 
       groups = id, type = "l", col = rgb(.1, .1, .1, .1),
       lwd = 2, layout = c(2, 1),
       ylab = "Neck pain (on visual analog scale)", xlab = "Examinations")
plot(p1)


## ----mtram-ordinalCont--------------------------------------------------------
library("ordinalCont")

## ----mtram-neck_ocm, cache = FALSE, results = "hide"--------------------------
neck_ocm <- ocm(vas ~ laser * time + (1 | id), data = pain_df, 
                scale = c(0, 1))


## ----mtram-neck_ocm_summary---------------------------------------------------
summary(neck_ocm)


## ----mtram-neck_ocm_ci--------------------------------------------------------
exp(cbind(coef(neck_ocm)[2:6], confint(neck_ocm)[2:6,]))


## ----tramME-neck--------------------------------------------------------------
library("tramME")
neck_ColrME <- ColrME(vas ~ laser * time + (1 | id), data = pain_df, 
                      bounds = c(0, 1), support = c(0, 1))


## ----tramME-neck_ci-----------------------------------------------------------
exp(coef(neck_ColrME))


## ----mtram-neck_Colr, cache = FALSE-------------------------------------------
neck_Colr <- Colr(vas ~ laser * time, data = pain_df, 
                  bounds = c(0, 1), support = c(0, 1),
                  extrapolate = TRUE)
neck_Colrmer <- mtram(neck_Colr, ~ (1 | id), data = pain_df)


## ----mtram-neck_Colr_distr, echo = FALSE, fig.height = 4, fig.width = 7-------
nd <- expand.grid(laser = unique(pain_df$laser),
                  time = unique(pain_df$time))
tmp <- as.mlt(neck_Colr)
coef(tmp)[] <- coef(neck_Colrmer)[1:12]
q <- 1:99/100
nd2 <- expand.grid(vas = q, laser = unique(pain_df$laser),
                   time = unique(pain_df$time))
# tmp <- as.mlt(neck_Colr) 
sd <- sqrt(coef(neck_Colrmer)[13]^2 + 1)
prb <- predict(tmp, newdata = nd, type = "distribution", q = q)
nd2$prob <- c(pnorm(qnorm(prb) / sd))
p2 <- xyplot(prob ~ vas | time, data = nd2, groups = laser, type = "l", 
             col = col, 
             layout = c(3, 1),
             xlab = "Neck pain (on visual analog scale)", 
             ylab = "Marginal distribution", 
             par.settings = simpleTheme(col=col),
             auto.key = list(text = levels(nd2$laser), 
                             points = FALSE, lines = TRUE, space = "top"))
plot(p2)

## M1
# neck_Colrmer <- mtram(neck_Colr, ~ (1 | id), data = pain_df, 
#                       Hessian = TRUE, standardise = FALSE)
# logLik(neck_Colrmer)
# 
# nd <- expand.grid(laser = unique(pain_df$laser),
#                   time = unique(pain_df$time))
# q <- 1:99/100
# nd2 <- expand.grid(vas = q, laser = unique(pain_df$laser),
#                    time = unique(pain_df$time))
# tmp <- as.mlt(neck_Colr)
# coef(tmp)[] <- coef(neck_Colrmer)[1:12]
# sd <- sqrt(coef(neck_Colrmer)[13]^2 + 1)
# prb <- predict(tmp, newdata = nd, type = "distribution", q = q)
# nd2$prob <- c(pnorm(qnorm(prb) / sd))
# p2 <- xyplot(prob ~ vas | time, data = nd2, groups = laser, type = "l", 
#              col = col, ylim = c(-0.05, 1.05),
#              layout = c(3, 1),
#              xlab = "Neck pain (on visual analog scale)", 
#              ylab = "Marginal distribution", 
#              par.settings = simpleTheme(col=col),
#              auto.key = list(text = levels(nd2$laser), 
#                              points = FALSE, lines = TRUE, space = "top"))
# plot(p2)


## ----mtram-neck_Colr-CI, echo = TRUE, eval=TRUE-------------------------------
S <- vcov(neck_Colrmer)
rbeta <- rmvnorm(10000, mean = coef(neck_Colrmer), sigma = S)
s <- rbeta[, ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
t(apply(rbeta[, 8:12], 2, function(x) {
  quantile(exp(x),prob = c(.025, .5, .975))}))


## ----mtram-neck_Colr-std_beta-------------------------------------------------
beta <- coef(neck_Colrmer)[8:12]
alpha <- coef(neck_Colrmer)[13]
(std_beta <- cbind(beta / sqrt(1 + alpha^2)))


## ----mtram-neck_Colr-ctr_mat--------------------------------------------------
ctr_mat <- matrix(c(1, 0, 0, 0, 0,
                    1, 0, 0, 1, 0,
                    1, 0, 0, 0, 1), nrow = 3, byrow = TRUE)
ctr_mat %*% std_beta


## ----mtram-neck_PImanual, eval=TRUE-------------------------------------------
(ci_emp <- t(apply(ctr_mat %*% t(rbeta[, 8:12]), 1, function(x) {
  quantile(x, prob = c(.025, .5, .975))})))

PI(-ci_emp, link = "logistic")


## ----neck_Colr-PI, echo = TRUE------------------------------------------------
nd <- expand.grid(time = unique(pain_df$time),
                  laser = unique(pain_df$laser))
neck_Colr_marg <- neck_Colr
neck_Colr_marg$coef <- coef(neck_Colrmer)[1:12] / 
                       sqrt(coef(neck_Colrmer)[13]^2 + 1)
(neck_Colr_PI <- PI(neck_Colr_marg, newdata = nd[1:3, ], 
                    reference = nd[4:6, ],
                    one2one = TRUE, conf.level = .95))[1:3, 1:3]


## ----mtram-CAO, echo = FALSE--------------------------------------------------
dir <- system.file("rda", package = "TH.data")
load(file.path(dir, "Primary_endpoint_data.rda"))


## ----mtram-CAO-plot, cache = FALSE, echo = FALSE------------------------------
ra <- sort(unique(CAOsurv$randarm))
st <- sort(unique(CAOsurv$strat_t))
sn <- sort(unique(CAOsurv$strat_n))
su <- c(1, 1700)
add <- c(0,  max(CAOsurv$iDFS[, "time2"]) - su[2])
ylim <- c(-.05, 1.05)
tmp <- as.mlt(Coxph(iDFS | 0 + strat_n:strat_t:randarm ~ 1, data = CAOsurv, 
                    support = su, add = add, log_first = TRUE))
nd <- expand.grid(strat_n = sn, strat_t = st, randarm = ra)
q <- mkgrid(tmp, 100)[[1]]
surv <- predict(tmp, newdata = nd, type = "survivor", q = q)
nd <- nd[rep(1:nrow(nd), each = nrow(surv)),]
nd$time <- q
nd$surv <- c(surv)
xyplot(surv ~ time | strat_t + strat_n, data = nd, groups = randarm, 
       type = "l", ylim = c(0, 1), col = col, ylab = "Probability",
       xlab = "Time (in days)",
       par.settings = simpleTheme(col=col),
       auto.key = list(text = levels(nd$randarm), 
                       points = FALSE, lines = TRUE, space = "top"))


## ----mtram-CAO_DFS------------------------------------------------------------
### convert "exact" event dates to interval-censoring (+/- one day)
tmp <- CAOsurv$iDFS
exact <- tmp[,3] == 1 
tmp[exact,2] <- tmp[exact,1] + 2
tmp[exact,1] <- pmax(tmp[exact,1] - 2, 0)
tmp[exact,3] <- 3
CAOsurv$iDFS2 <- tmp


## ----mtram-CAO_SR, cache = TRUE-----------------------------------------------
CAO_SR <- Survreg(iDFS2 ~ randarm, data = CAOsurv)
CAO_SR_mtram <- mtram(CAO_SR, ~ (1 | Block), data = CAOsurv)
logLik(CAO_SR_mtram)
(cf <- coef(CAO_SR_mtram))
(OR <- exp(-cf["randarm5-FU + Oxaliplatin"] / sqrt(cf["gamma1"]^2 + 1)))


## ----mtram-CAO-CI-------------------------------------------------------------
S <- vcov(CAO_SR_mtram)
# sqrt(diag(S))
rbeta <- rmvnorm(10000, mean = coef(CAO_SR_mtram), 
                 sigma = S)
s <- rbeta[, ncol(rbeta)]
rbeta <- rbeta[, -ncol(rbeta)] / sqrt(s^2 + 1)
quantile(exp(-rbeta[, ncol(rbeta)]), prob = c(.025, .5, .975))


## ----mtram-CAO_SR_2, cache = TRUE---------------------------------------------
CAO_SR_2 <- Survreg(iDFS2 | 0 + strat_n:strat_t ~ randarm, data = CAOsurv)
CAO_SR_2_mtram <- mtram(CAO_SR_2, ~ (1 | Block), data = CAOsurv)
logLik(CAO_SR_2_mtram)
(cf <- coef(CAO_SR_2_mtram))
(OR_2 <- exp(-cf["randarm5-FU + Oxaliplatin"] / sqrt(cf["gamma1"]^2 + 1)))


## ----mtram-CAO-CI-2, echo = FALSE---------------------------------------------
S <- vcov(CAO_SR_2_mtram)
rbeta <- rmvnorm(10000, mean = coef(CAO_SR_2_mtram), 
                 sigma = S)
s <- rbeta[, ncol(rbeta)]
rbeta <- rbeta[, -ncol(rbeta)] / sqrt(s^2 + 1)
quantile(exp(-rbeta[, ncol(rbeta)]), prob = c(.025, .5, .975))


## ----mtram-CAO_Cox_2, cache = TRUE--------------------------------------------
CAO_Cox_2 <- Coxph(iDFS2 | 0 + strat_n:strat_t ~ randarm, data = CAOsurv, 
                   support = c(1, 1700), log_first = TRUE, order = 4)
logLik(CAO_Cox_2)
CAO_Cox_2_mtram <- mtram(CAO_Cox_2, ~ (1 | Block), data = CAOsurv)
logLik(CAO_Cox_2_mtram)
coef(CAO_Cox_2_mtram)


## ----mtram-CAO-CI-3, echo = FALSE---------------------------------------------
S <- vcov(CAO_Cox_2_mtram)
rbeta <- rmvnorm(10000, mean = coef(CAO_Cox_2_mtram), 
                 sigma = S)
s <- rbeta[,ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
quantile(exp(rbeta[, ncol(rbeta)]), prob = c(.025, .5, .975))


## ----mtram-CAO_PI-------------------------------------------------------------
nd <- CAOsurv[1:2, ]
tmp <- CAO_Cox_2
tmp$coef <- coef(CAO_Cox_2_mtram)[-22] / sqrt(coef(CAO_Cox_2_mtram)[22]^2 + 1)
(CAO_Cox_PI <- PI(tmp, newdata = nd[2, ], reference = nd[1, ],
                  one2one = TRUE, conf.level = .95))[1, ]


## ----mtram-CAO_PI_man---------------------------------------------------------
ci_man <- quantile(-rbeta[, ncol(rbeta)], prob = c(.025, .5, .975))
(CAO_Cox_PIm <- PI(ci_man, link = "minimum extreme value"))


## ----tramME-CAO_SR, cache = TRUE----------------------------------------------
CAO_Cox_2_tramME <- CoxphME(iDFS2 | 0 + strat_n:strat_t ~ randarm + (1 | Block), 
                            data = CAOsurv, log_first = TRUE)


## ----tramME-CAO_SR-hr, cache = TRUE-------------------------------------------
exp(coef(CAO_Cox_2_tramME))
exp(confint(CAO_Cox_2_tramME, parm = "randarm5-FU + Oxaliplatin", 
            estimate = TRUE))


## ----echo=FALSE, eval=FALSE---------------------------------------------------
# sqrt(VarCorr(CAO_Cox_2_tramME)$Block$var)
# coef(CAO_Cox_2_mtram)["gamma1"]


## ----mtram-CAO-coxme, echo = FALSE, eval = FALSE------------------------------
# library("coxme")
# m <- coxme(DFS ~ randarm + (1 | Block), data = CAOsurv)
# summary(m)
# sd <- sqrt(diag(vcov(m)))
# exp(coef(m) + c(-1, 0, 1) * qnorm(.975) * sd)


## ----sim, eval = FALSE--------------------------------------------------------
# source(system.file("simulations", "mtram_sim.R", package = "tram"), echo = TRUE)


## ----mtram-sessionInfo, echo = FALSE, results = "hide"------------------------
sessionInfo()


## ----mtram-funs, echo = FALSE, results = "hide"-------------------------------
if (file.exists("packages.bib")) file.remove("packages.bib")
pkgversion <- function(pkg) {
    pkgbib(pkg)
    packageDescription(pkg)$Version
}
pkgbib <- function(pkg) {
    x <- citation(package = pkg, auto = TRUE)[[1]]
    b <- toBibtex(x)
    b <- gsub("Buehlmann", "B{\\\\\"u}hlmann", b)
    b[1] <- paste("@Manual{pkg:", pkg, ",", sep = "")
    if (is.na(b["url"])) {
        b[length(b)] <- paste("   URL = {http://CRAN.R-project.org/package=",
                              pkg, "}", sep = "")
        b <- c(b, "}")
    }
    cat(b, sep = "\n", file = "packages.bib", append = TRUE)
}

pkg <- function(pkg) {
    vrs <- try(pkgversion(pkg))
    if (inherits(vrs, "try-error")) return(NA)
    paste("\\\\pkg{", pkg, "} \\\\citep[version~",
          vrs, ",][]{pkg:", pkg, "}", sep = "")
}

pkg("mlt")
pkg("tram")
pkg("SparseGrid")
cat(c("@Manual{vign:mlt.docreg,",
             "    title = {Most Likely Transformations: The mlt Package},",
             "    author = {Torsten Hothorn},",
             paste("    year = ", substr(packageDescription("mlt.docreg")$Date, 1, 4), ",", sep = ""),
             paste("    note = {R package vignette version ", packageDescription("mlt.docreg")$Version, "},", sep = ""),
             "    url = {https://CRAN.R-project.org/package=mlt.docreg},",
             "}"), file = "packages.bib", append = TRUE, sep = "\n")

