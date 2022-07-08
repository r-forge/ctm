# Demo cotram

# Windows diffs...
options(digits = 3)

set.seed(24101968)
library(tramvs)
library(cotram)

data("birds", package = "TH.data")
birds$noise <- rnorm(nrow(birds), sd = 10)

# Estimate support sice via HBIC
res <- tramvs(SG5 ~ AOT + AFS + GST + DBH + DWC + LOG + noise, data = birds,
              modFUN = cotram)
plot(res, type = "b")
plot(res, which = "path")

# Active set
support(res)
coef(res, best_only = TRUE)
coef(res, best_only = FALSE, with_baseline = TRUE)
coef(res, best_only = TRUE, with_baseline = TRUE)
