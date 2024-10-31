
library("tram")
require("lme4")

## ----mtram-sleep_lmer, cache = FALSE------------------------------------------
sleep_lmer <- lmer(Reaction ~ Days + (Days | Subject), 
                   data = sleepstudy, REML = FALSE)


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
