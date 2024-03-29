
library("tram")
library("survival")
library("trtf")

### Windows diffs...
options(digits = 3)

chk <- function(x, y, tol = 1e-1) 
    stopifnot(isTRUE(all.equal(x, y, check.attributes = FALSE, tol = tol)))

data("GBSG2", package = "TH.data")

cmod <- coxph(Surv(time, cens) ~ progrec + pnodes + strata(horTh, tgrade),
               data = GBSG2)
Cmod <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ progrec + pnodes, 
              data = GBSG2)

chk(coef(cmod), coef(Cmod))
chk(diag(vcov(cmod)), diag(vcov(Cmod)))

Cmod_lf <- Coxph(Surv(time, cens) | 0 + horTh:tgrade ~ progrec + pnodes, 
                data = GBSG2, log_first = TRUE)

chk(coef(cmod), coef(Cmod_lf))
chk(diag(vcov(cmod)), diag(vcov(Cmod_lf)))

cmod_2 <- coxph(Surv(time, cens) ~ ., data = GBSG2)
Cmod_2 <- Coxph(Surv(time, cens) ~ ., data = GBSG2)

chk(coef(cmod_2), coef(Cmod_2))
chk(diag(vcov(cmod_2)), diag(vcov(Cmod_2)))

cmod <- Coxph(Surv(time, cens) ~ horTh, data = GBSG2)
(tmod <- trafotree(cmod, formula = Surv(time, cens) ~ horTh | ., data = GBSG2))
logLik(tmod)

### check residuals
library("sandwich")
GBSG2$y <- with(GBSG2, Surv(time, cens))
ORDER <- 12
### indirect computation: score wrt to int
GBSG2$int <- 1
m <- Coxph(y ~ int, data = GBSG2, fixed = c("int" = 0), 
           LRtest = FALSE, order = ORDER)
m1 <- mlt(m$model, data = GBSG2, dofit = FALSE)
coef(m1) <- coef(as.mlt(m))
LR1 <- estfun(m1)[, length(coef(m1))]

### direct computation
m2 <- Coxph(y ~ 1, data = GBSG2, 
            LRtest = FALSE, order = ORDER)
LR2 <- resid(as.mlt(m2))

chk(LR1, LR2)

## interval-censoring
load(system.file("rda", "Primary_endpoint_data.rda", package = "TH.data"))

## [tram] Parametric distribution-free PH model (interval censoring)
mci <- Coxph(iDFS ~ randarm, data = CAOsurv, log_first = TRUE)
logLik(mci)
