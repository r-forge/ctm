set.seed(24101968)

# Windows diffs...
options(digits = 3)
library(tramvs)

N <- 1e2
P <- 5
nz <- 3
beta <- rep(c(1, 0), c(nz, P - nz))
X <- matrix(rnorm(N * P), nrow = N, ncol = P)
Y <- 1 + X %*% beta + rnorm(N)

dat <- data.frame(y = Y, x = X)

# Mandatory noise covariate
tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
              mandatory = y ~ x.5)

# Mandatory noise covariate in shift and scale
tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
              mandatory = y ~ x.5 | x.5)

# Mandatory noise covariate in scale only
tramvs(y ~ . | x.5, data = dat, modFUN = Lm,
              mandatory = y ~ 1 | x.5)
