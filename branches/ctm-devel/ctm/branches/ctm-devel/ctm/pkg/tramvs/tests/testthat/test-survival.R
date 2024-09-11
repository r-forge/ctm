### Test survival

test_that("survival", {
  set.seed(24101968)
  library("tramvs")

  if (require("TH.data") & require("abess") & require("survival")) {
    data("GBSG2", package = "TH.data")
    GBSG2$surv <- Surv(GBSG2$time, GBSG2$cens)

    res <- tramvs(surv ~ horTh + age + menostat + tsize + tgrade +
                    pnodes + progrec + estrec, data = GBSG2, modFUN = Coxph)
    res_abess <- abess(model.matrix(~ horTh + age + menostat + tsize + tgrade +
                                      pnodes + progrec + estrec, data = GBSG2)[, -1],
                       y = with(GBSG2, Surv(time, cens)),
                       data = GBSG2, family = "cox")

    # M1 diffs
    round(coef(res), 2)
    round(coef(res_abess)[,-1], 2)

    # Active set
    expect_equal(support(res), c("tgrade.L", "pnodes", "progrec"))
  }
})
