# Demo continuous / sparse

# Windows diffs...
options(digits = 3)

set.seed(24101968)
library(tramvs)
library(abess)

N <- 1e2
P <- 10
nz <- 3
beta <- rep(c(3, 0), c(nz, P - nz))
X <- matrix(rnorm(N * P), nrow = N, ncol = P)
Y <- 1 + X %*% beta + rnorm(N)

dat <- data.frame(y = Y, x = X)
res <- tramvs(y ~ ., data = dat, modFUN = Lm)
res_abess <- abess(y ~ ., data = dat, family = "gaussian")

as(as.matrix(coef(res, as.lm = TRUE)), "sparseMatrix")
coef(res_abess)[,-1]

# S3 methods

print(res)
summary(res)
plot(res, type = "b")
plot(res, which = "path")
logLik(res)
SIC(res)
coef(res)
predict(res, which = "distribution", type = "trafo")
simulate(res)
residuals(res)

# Active set
support(res)
