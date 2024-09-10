### Test Weibull

test_that("weibull", {
  set.seed(24101968)
  library("tramvs")

  if (require("TH.data") & require("survival")) {
    data("GBSG2", package = "TH.data")

    ## using the alias SurvregVS
    GBSG2$surv <- with(GBSG2, Surv(time, cens))
    res <- SurvregVS(surv ~ horTh + age + menostat + tsize + tgrade +
                       pnodes + progrec + estrec, data = GBSG2)

    coef(res)

    # Active set
    expect_equal(support(res), c("tgrade.L", "pnodes", "progrec"))
  }
})
