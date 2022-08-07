# Demo count / sparse

# Windows diffs...
old <- options(digits = 3)

set.seed(24101968)
library(tramvs)
library(cotram)

N <- 1e2
P <- 5
nz <- 3
beta <- rep(c(1, 0), c(nz, P - nz))
X <- matrix(abs(rnorm(N * P)), nrow = N, ncol = P)
Y <- as.integer(1 + X %*% beta + abs(rnorm(N)))

dat <- data.frame(y = Y, x = X)
res <- tramvs(y ~ ., data = dat, modFUN = cotram)

# Active set
support(res)

options(old)
