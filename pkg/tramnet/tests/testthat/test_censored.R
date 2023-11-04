# Tests for COLR models

test_that("Colr and Coxph", {
  set.seed(0)

  ## Depencies
  library("sandwich")

  if (require("survival") & require("TH.data")) {
    ## Exact and Right censored
    data("GBSG2", package = "TH.data")
    GBSG2$surv <- with(GBSG2, Surv(time, cens))

    modCOLR <- ColrNET(surv ~ horTh, data = GBSG2, tram_args = list(
      log_first = TRUE, order = 4, support = c(1e-12, max(GBSG2$time))))
    yCOLRb <- Colr(surv ~ horTh, data = GBSG2, log_first = TRUE, order = 4,
                   support = c(1e-12, max(GBSG2$time)))
    expect_lt(max(abs(coef(yCOLRb, with_baseline = FALSE) -
                        coef(modCOLR, with_baseline = FALSE) /
                        sd(as.numeric(GBSG2$horTh)))),
              0.01)
    expect_lt(abs(logLik(modCOLR) + modCOLR$result$value), 0.1)
    expect_no_error({
      ## methods
      logLik(modCOLR, newdata = GBSG2[2, ])
      coef(modCOLR, tol = 0, with_baseline = TRUE)
      logLik(modCOLR)
      c(resid(modCOLR)[1:10])
      predict(modCOLR, type = "distribution", q = 1)[, 1:10]
      predict(modCOLR, type = "quantile", prob = 0.5)[, 1:10]
      unclass(simulate(modCOLR)[1:5,])
      as.data.frame(head(estfun(modCOLR)))
      plot(modCOLR, type = "survivor")
      plot(modCOLR, type = "density", K = 120)
      print(modCOLR)
    })
  }
})
