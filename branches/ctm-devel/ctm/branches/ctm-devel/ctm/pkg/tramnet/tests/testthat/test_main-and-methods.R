# Tests for Box-Cox type regularized regression models

test_that("Regularized Box-Cox models", {
  ## Set up and model fit
  data("cars", package = "datasets")
  cars <- as.data.frame(scale(cars))
  mt <- BoxCoxNET(speed ~ dist, data = cars, alpha = 0, lambda = 0)
  m2 <- BoxCox(speed ~ dist, data = cars)

  expect_lt(max(abs(coef(mt, with_baseline = FALSE) -
                      coef(m2, with_baseline = FALSE))), 1e-5)

  expect_equal(logLik(mt)[1], logLik(m2)[1])

  expect_no_error(cvl_tramnet(mt))

  expect_no_error({
    ## methods
    logLik(mt, newdata = cars[2, ])
    coef(mt, tol = 0, with_baseline = TRUE)
    c(residuals(mt)[1:10])
    predict(mt, type = "distribution", q = 1)[, 1:10]
    as.double(predict(mt, type = "quantile", prob = 0.5))
    as.double(simulate(mt)[1:5,])
    as.data.frame(head(estfun(mt)))
    plot(mt, type = "survivor")
    plot(mt, type = "density", K = 120)
    print(mt)
  })
})

test_that("Constraints", {
  ## Test for additional inequality constraints on beta
  data("cars", package = "datasets")
  m2 <- BoxCox(speed ~ dist, data = cars, constraints = c("dist <= 0"))
  lhs <- attr(model.matrix(m2), "constraint")$ui
  rhs <- attr(model.matrix(m2), "constraint")$ci
  mt <- BoxCoxNET(speed ~ dist, data = cars, alpha = 0, lambda = 0, constraints = list(lhs, rhs))

  expect_lt(
    max(abs(coef(mt, with_baseline = FALSE) -
              coef(m2, with_baseline = FALSE)[-2])),
    1e-5
  )

  expect_equal(logLik(mt)[1], logLik(m2)[1], tolerance = 1e-3)
})

test_that("Alias", {
  expect_no_error({
    data("cars", package = "datasets")
    LmNET(speed ~ dist, data = cars, alpha = 0, lambda = 0)
    SurvregNET(speed ~ dist, data = cars, alpha = 0, lambda = 0)
    LehmannNET(speed ~ dist, data = cars, alpha = 0, lambda = 0)
    CoxphNET(speed ~ dist, data = cars, alpha = 0, lambda = 0)
  })
})

test_that("Stratified", {
  dat <- data.frame(y = runif(100), s = factor(rep(c(1, 2), each = 50)))
  x <- scale(matrix(rnorm(100 * 20, mean = 0, sd = 1), nrow = 100))
  colnames(x) <- paste0("X", 1:20)
  y2 <- Lm(y | 0 + s ~ 1, data = dat)
  expect_no_error(tramnet(y2, x, lambda = 8, alpha = 1))
})
