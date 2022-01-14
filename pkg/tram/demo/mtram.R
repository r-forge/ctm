## ----mtram-setup, echo = FALSE, results = "hide", message = FALSE, warning = FALSE----
set.seed(290875)

pkgs <- sapply(c("mlt", "survival", "tram", "lme4", "gridExtra", 
         "lattice", "latticeExtra", "mvtnorm", "ordinalCont"), require, char = TRUE)

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
library("colorspace")
col <- diverge_hcl(2, h = c(246, 40), c = 120, l = c(65, 90), alpha = .75)


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


## ----mtram-vignette, eval = FALSE---------------------------------------------
## install.packages("tram")
## vignette("mtram", package = "tram")


## ----mtram-sleep-plot, echo = FALSE-------------------------------------------
xyplot(Reaction ~ Days | Subject, data = sleepstudy, 
       xlab = "Days of sleep deprivation", ylab = "Average reaction time (in ms)")


## ----mtram-sleep_lmer, cache = TRUE-------------------------------------------
library("lme4")
sleep_lmer <- lmer(Reaction ~ Days + (Days | Subject), 
                   data = sleepstudy, REML = FALSE)


## ----mtram-sleep_mtram, cache = TRUE------------------------------------------
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


## ----mtram-sleep-interval, cache = TRUE---------------------------------------
sleep_LM_I <- Lm(Reaction_I ~ Days, data = sleepstudy)
sleep_LMmer_I <- mtram(sleep_LM_I, ~ (Days | Subject), data = sleepstudy)


## ----mtram-sleep-interval-results---------------------------------------------
logLik(sleep_LMmer_I)
coef(sleep_LMmer_I)
coef(sleep_LMmer)


## ----mtram-sleep_BoxCox, cache = TRUE-----------------------------------------
sleep_BC <- BoxCox(Reaction ~ Days, data = sleepstudy)
sleep_BCmer <- mtram(sleep_BC, ~ (Days | Subject), data = sleepstudy, 
                     Hessian = TRUE)
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
meth <- c("Normal linear mixed-effects model", "Non-norma linear transformation model")
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
                  points = FALSE, lines = TRUE))


## ----mtram-sleep_corr---------------------------------------------------------
cov2cor(sleep_BCmer$G)


## ----mtram-sleep_vcov---------------------------------------------------------
VC <- solve(sleep_BCmer$Hessian)
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


## ----mtram-toenail-plot, echo = FALSE-----------------------------------------
data("toenail", package = "HSAUR3")
rlev <- levels(toenail$patientID)[xtabs(~ patientID, 
                                        data = toenail) == 1]
toenail <- subset(toenail, !patientID %in% rlev)
toenail$patientID <- toenail$patientID[, drop = TRUE]
layout(matrix(1:2, ncol = 2))
trt <- levels(toenail$treatment)
cdplot(outcome ~ time, data = subset(toenail, treatment == trt[1]),
       main = trt[1], xlab = "Time", ylab = "Toe nail infection")
cdplot(outcome ~ time, data = subset(toenail, treatment == trt[2]),
       main = trt[2], xlab = "Time", ylab = "")


## ----mtram-toenail_glmer_RI, cache = TRUE-------------------------------------
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


## ----mtram-toenail_mtram_RI, cache = TRUE-------------------------------------
m <- ctm(as.basis(~ outcome, data = toenail), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Normal")
toenail_probit <- mlt(m, data = toenail, 
                      fixed = c("outcomemoderate or severe" = 0))
toenail_mtram_RI <- 
    mtram(toenail_probit, ~ (1 | patientID), 
          data = toenail, Hessian = TRUE)
logLik(toenail_mtram_RI)
coef(toenail_mtram_RI)


## ----mtram-toenail-hessian----------------------------------------------------
vcov(toenail_glmer_RI_2)
solve(toenail_mtram_RI$Hessian)[1:4, 1:4]


## ----mtram-toenail_glmer_RS, cache = TRUE-------------------------------------
toenail_glmer_RS <- 
    glmer(outcome ~ treatment * time + (1 + time | patientID),
          data = toenail, family = binomial(link = "probit"))
summary(toenail_glmer_RS)
toenail_glmer_RS@theta


## ----mtram-toenail_mtram_RS, cache = TRUE-------------------------------------
toenail_mtram_RS <- 
    mtram(toenail_probit, ~ (1 + time | patientID), 
          data = toenail)
logLik(toenail_mtram_RS)
coef(toenail_mtram_RS)


## ----mtram-toenail_logit, cache = TRUE----------------------------------------
m <- ctm(as.basis(~ outcome, data = toenail), 
         shifting = ~ treatment * time, 
         data = toenail, todistr = "Logistic")
toenail_logit <- mlt(m, data = toenail, 
                     fixed = c("outcomemoderate or severe" = 0))
toenail_mtram_logit <- mtram(toenail_logit, ~ (1 | patientID), 
                              data = toenail)
toenail_mtram_logit_s <- mtram(toenail_logit, ~ (1 | patientID), 
                               data = toenail, standardise = TRUE, 
                               Hessian = TRUE)


## ----mtram-toenail_marginal_logit_s-------------------------------------------
tmp <- toenail_logit
cf <- coef(tmp)
cf <- cf[names(cf) != "outcomemoderate or severe"]
sdrf <- rev(coef(toenail_mtram_logit_s))[1]
cf <- coef(toenail_mtram_logit_s)[names(cf)] / sqrt(sdrf^2 + 1)
cf <- c(cf[1], "outcomemoderate or severe" = 0, cf[-1])
coef(tmp) <- cf
time <- 0:180/10
treatment <- sort(unique(toenail$treatment))
nd <- expand.grid(time = time, treatment = treatment)
nd$prob_logit_s <- predict(tmp, newdata = nd, type = "distribution")[1,]
nd$odds <- exp(predict(tmp, newdata = nd, type = "trafo")[1,])


## ----mtram-toenail_OR_2, dev = "png", cache = TRUE, echo = FALSE, dpi = 300----
X <- model.matrix(~ treatment * time, data = nd)
rbeta <- rmvnorm(10000, mean = coef(toenail_mtram_logit_s), 
                 sigma = solve(toenail_mtram_logit_s$Hessian))
s <- rbeta[,ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
odds <- exp(X %*% t(rbeta))
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
nd2 <- nd[rep(1:nrow(nd), 3),]
nd2$prob <- c(nd$prob_probit, nd$prob_logit_s, nd$prob_logit)
lev <- c("Probit (M1) = (M2)", "Logit (M1)", "Logit (M2)")
nd2$model <- rep(factor(lev, labels = lev, levels = lev), each = nrow(nd))

xyplot(prob ~ time | model, data = nd2, group = treatment, ylim = c(0, 1), 
       xlab = "Time", 
       par.settings = simpleTheme(col=col),
       auto.key = list(text = levels(nd2$treatment), 
                       points = FALSE, lines = TRUE), 
       col = col, type = "l", ylab = "Probability (none or mild)")


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


## ----mtram-neck_ocm, cache = TRUE, results = "hide"---------------------------
library("ordinalCont")
neck_ocm <- ocm(vas ~ laser * time + (1 | id), data = pain_df, 
                scale = c(0, 1))


## ----mtram-neck_ocm_summary---------------------------------------------------
summary(neck_ocm)


## ----mtram-neck_ocm_ci--------------------------------------------------------
exp(cbind(coef(neck_ocm)[2:6], confint(neck_ocm)[2:6,]))


## ----mtram-neck_Colr, cache = TRUE--------------------------------------------
neck_Colr <- Colr(vas ~ laser * time, data = pain_df, 
                  bounds = c(0, 1), support = c(0, 1),
                  extrapolate = TRUE)
neck_Colrmer <- mtram(neck_Colr, ~ (1 | id), data = pain_df, 
                      Hessian = TRUE)
logLik(neck_Colrmer)


## ----mtram-neck_Colr_distr, echo = FALSE, fig.height = 4, fig.width = 7-------
nd <- expand.grid(laser = unique(pain_df$laser),
                  time = unique(pain_df$time))
tmp <- as.mlt(neck_Colr)
coef(tmp)[] <- coef(neck_Colrmer)[1:12]
q <- 1:99/100
nd2 <- expand.grid(vas = q, laser = unique(pain_df$laser),
                   time = unique(pain_df$time))
tmp <- as.mlt(neck_Colr)
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
                       points = FALSE, lines = TRUE))
plot(p2)


## ----mtram-CAO, echo = FALSE--------------------------------------------------
dir <- system.file("rda", package = "TH.data")
load(file.path(dir, "Primary_endpoint_data.rda"))


## ----mtram-CAO-plot, cache = TRUE, echo = FALSE-------------------------------
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
                       points = FALSE, lines = TRUE))


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
CAO_SR_mtram_s <- mtram(CAO_SR, ~ (1 | Block), data = CAOsurv,
                        standardise = TRUE, Hessian = TRUE)
logLik(CAO_SR_mtram_s)
(cf <- coef(CAO_SR_mtram_s))
(OR <- exp(-cf["randarm5-FU + Oxaliplatin"] / sqrt(cf["gamma1"]^2 + 1)))


## ----mtram-CAO-CI-------------------------------------------------------------
S <- solve(CAO_SR_mtram_s$Hessian)
sqrt(diag(S))
rbeta <- rmvnorm(10000, mean = coef(CAO_SR_mtram_s), 
                 sigma = S)
s <- rbeta[,ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
quantile(exp(-rbeta[, ncol(rbeta)]), prob = c(.025, .5, .975))


## ----mtram-CAO_SR_2, cache = TRUE---------------------------------------------
CAO_SR_2 <- Survreg(iDFS2 | 0 + strat_n:strat_t ~ randarm, data = CAOsurv)
CAO_SR_2_mtram_s <- mtram(CAO_SR_2, ~ (1 | Block), data = CAOsurv,
                          standardise = TRUE, Hessian  = TRUE)
logLik(CAO_SR_2_mtram_s)
(cf <- coef(CAO_SR_2_mtram_s))
(OR_2 <- exp(-cf["randarm5-FU + Oxaliplatin"] / sqrt(cf["gamma1"]^2 + 1)))


## ----mtram-CAO-CI-2, echo = FALSE---------------------------------------------
H <- CAO_SR_2_mtram_s$Hessian
S <- solve(H + .01 * diag(nrow(H)))
sqrt(diag(S))
rbeta <- rmvnorm(10000, mean = coef(CAO_SR_2_mtram_s), 
                 sigma = S)
s <- rbeta[,ncol(rbeta)]
rbeta <- rbeta[,-ncol(rbeta)] / sqrt(s^2 + 1)
quantile(exp(-rbeta[, ncol(rbeta)]), prob = c(.025, .5, .975))


## ----mtram-CAO-coxme----------------------------------------------------------
library("coxme")
m <- coxme(DFS ~ randarm + (1 | Block), data = CAOsurv)
summary(m)
sd <- sqrt(diag(vcov(m)))
exp(coef(m) + c(-1, 0, 1) * qnorm(.975) * sd)


## ----mtram-CHFLS, echo = FALSE------------------------------------------------
library("TH.data")
load(file.path(path.package(package="TH.data"), "rda", "CHFLS.rda"))

### choose neccessary variables
org <- chfls1[, c("REGION6", "ZJ05", "ZJ06", "A35", "ZJ07", "ZJ16M", "INCRM",
                  "JK01", "JK02", "JK20", "HY04", "HY07", "A02", "AGEGAPM", 
                  "A07M", "A14", "A21", "A22M", "A23", "AX16", "INCAM", "SEXNOW", "ZW04")]

names(org) <- c("Region",
                "Rgender",               ### gender of respondent
                "Rage",                  ### age of respondent
                "RagestartA",            ### age of respondent at beginning of relationship with partner A
                "Redu",                  ### education of respondent
                "RincomeM",              ### rounded monthly income of respondent
                "RincomeComp",           ### inputed monthly income of respondent
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

olevels <- c("never", "rarely", "sometimes", "often", "always")
orgA <- subset(org, Rgender == "female" & Rhomosexual != "yes" & orgasm %in% olevels)

orgA$orgasm <- ordered(as.character(orgA$orgasm),
        levels = c("never", "rarely", "sometimes", "often", "always"))

orgA$Redu <- factor(as.character(orgA$Redu),
        levels = c("univ/grad", "j col", "up mid", "low mid", "primary", "no school"))
levels(orgA$Redu) <-  c("univ", "jcol", "upmid", "lowmid", "primary", "noschool")

orgA$Aedu <- factor(as.character(orgA$Aedu),
        levels = c("univ/grad", "j col", "up mid", "low mid", "primary", "no school"))

orgA$Rhappy <- factor(as.character(orgA$Rhappy),
        levels = c("v unhappy", "not too", "relatively", "very"))

orgA$Rhealth <- factor(as.character(orgA$Rhealth),
        levels = c("poor", "not good", "fair", "good", "excellent"))

orgA$Region <- factor(as.character(orgA$Region),
        levels = c("CentralW", "Northeast", "North", "InlandS", "CoastalE", "CoastalS"))

orgA$AincomeSD <- orgA$AincomeComp/sd(orgA$AincomeComp)
orgA$AheightSD <- orgA$Aheight/sd(orgA$Aheight)
orgA$RageSD <- orgA$Rage/sd(orgA$Rage)
orgA$edudiff <- as.numeric(orgA$Aedu) - as.numeric(orgA$Redu)
orgA$edudiffSD <- orgA$edudiff/sd(orgA$edudiff, na.rm=TRUE)
orgA$wealthdiff <- orgA$RincomeComp - orgA$AincomeComp
orgA$wealthdiffSD <- orgA$wealthdiff/sd(orgA$wealthdiff, na.rm=TRUE)
orgA$RAdurationSD <- orgA$RAduration/sd(orgA$RAduration, na.rm=TRUE)
orgAtmp <- orgA[, c("orgasm", "AincomeSD", "AheightSD", "RAdurationSD",
                 "RageSD", "edudiffSD", "wealthdiffSD", "Redu", "Rhealth",
                 "Rhappy", "Region")]
cc <- complete.cases(orgAtmp)
orgAcc <- subset(orgA, cc)


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
#pkg("lme4")
cat(c("@Manual{vign:mlt.docreg,",
             "    title = {Most Likely Transformations: The mlt Package},",
             "    author = {Torsten Hothorn},",
             paste("    year = ", substr(packageDescription("mlt.docreg")$Date, 1, 4), ",", sep = ""),
             paste("    note = {R package vignette version ", packageDescription("mlt.docreg")$Version, "},", sep = ""),
             "    url = {https://CRAN.R-project.org/package=mlt.docreg},",
             "}"), file = "packages.bib", append = TRUE, sep = "\n")

