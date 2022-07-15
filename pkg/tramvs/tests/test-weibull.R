# Demo Weibull

# Windows diffs...
old <- options(digits = 3)

set.seed(24101968)
library(tramvs)
library(abess)
library(survival)

data("GBSG2", package = "TH.data")

## using the alias SurvregVS
GBSG2$surv <- with(GBSG2, Surv(time, cens))
res <- SurvregVS(surv ~ horTh + age + menostat + tsize + tgrade +
                   pnodes + progrec + estrec, data = GBSG2)

coef(res)

# Active set
support(res)

options(old)
