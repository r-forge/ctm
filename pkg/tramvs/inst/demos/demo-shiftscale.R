# Shift-scale transformation models
# TH, SS, LK
# Jan 2022

set.seed(29)

# Deps --------------------------------------------------------------------

library(tramvs)

# Params ------------------------------------------------------------------

N <- 1000
p <- 4
z <- rnorm(N)
x <- matrix(runif(N * p), ncol = p)
y <- (z + x[,1] * 2) / sqrt(exp(x[,2]))
d <- data.frame(y = y, x = x)
### 0, 1, 2, 1
coef(as.mlt(Lm(y ~ x.1 | x.2, data = d)))

fm <- paste(colnames(d)[-1], collapse = "+")
fm <- as.formula(paste("y ~ ", fm, "|", fm))

m0 <- LmVS(fm, data = d)
m0

coef(m0)
