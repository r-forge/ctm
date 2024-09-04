### Mandatory sets
test_that("mandatory", {
  set.seed(24101969)
  library("tramvs")

  N <- 1e2
  P <- 5
  nz <- 3
  beta <- rep(c(1, 0), c(nz, P - nz))
  X <- matrix(rnorm(N * P), nrow = N, ncol = P)
  Y <- 1 + X %*% beta + rnorm(N)

  dat <- data.frame(y = Y, x = X)

  # Mandatory noise covariate
  expect_no_error({
    tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
           mandatory = y ~ x.5)
  })

  # Mandatory noise covariate in shift and scale
  expect_no_error({
    tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
           mandatory = y ~ x.5 | x.5)
  })

  # Mandatory noise covariate in scale only
  expect_no_error({
    tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
           mandatory = y ~ 1 | x.5)
  })
})
