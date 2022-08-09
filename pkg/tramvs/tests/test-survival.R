# Demo survival

# Windows diffs...
old <- options(digits = 3)

set.seed(24101968)
library(tramvs)
library(abess)
library(survival)

data("GBSG2", package = "TH.data")

res <- tramvs(Surv(time, cens) ~ . - time - cens, data = GBSG2, modFUN = Coxph)
res_abess <- abess(model.matrix(~ . - time - cens, data = GBSG2)[, -1],
                   y = with(GBSG2, Surv(time, cens)),
                   data = GBSG2, family = "cox")

# M1 diffs
round(coef(res), 2)
round(coef(res_abess)[,-1], 2)

# Active set
support(res)

options(old)
